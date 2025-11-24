"""zmq functions for FEMAG

"""
import re
import pathlib
import threading
import json
import time
import logging
try:
    import zmq
except ImportError:
    pass

numpat = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)\s*')
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARN)

class ProtFile:
    def __init__(self, dirname, num_cur_steps):
        self.size = 0
        self.looplen = 0
        self.cur_steps = [1, num_cur_steps]
        self.n = 0
        self.num_loops = 0
        import platform
        self.dirname = dirname
        self.name = 'samples'

    def percent(self):
        if self.looplen > 0:
            return min(100 * self.n / self.looplen, 100)
        return 0

    def update(self):
        if not self.dirname:
            return ''
        p = list(pathlib.Path(self.dirname).glob('*.PROT'))
        if p:
            buf = ''
            if self.size < p[0].stat().st_size:
                with p[0].open() as fp:
                    fp.seek(self.size)
                    buf = fp.read()
            return self.append(buf)
        return ''

    def append(self, buf):
        self.size += len(buf)
        for line in [l.strip() for l in buf.split('\n') if l]:
            if line.startswith('Loop'):
                self.n = 0
                try:
                    cur_steps = self.cur_steps[self.num_loops]
                except IndexError:
                    cur_steps = 1
                x0, x1, dx, nbeta = [float(f)
                                     for f in re.findall(numpat, line)][:4]
                move_steps = round((x1-x0)/dx+1)
                beta_steps = int(nbeta)
                self.looplen = cur_steps*beta_steps*move_steps
                self.num_loops += 1
            elif (line.startswith('Cur') or
                  line.startswith('Id')):
                self.n += 1
            elif line.startswith('Number movesteps Fe-Losses'):
                return f'{self.percent():3.1f}%' # 100%
            elif line.startswith('begin'):
                self.name = line.split()[1].strip()

        return f'{self.percent():3.1f}%'  # {self.n}/{self.looplen}'

class SubscriberTask(threading.Thread):
    '''Subscribes data from all femag subprocesses

       topics with
       progress, xydata and all entries of "simple_topic" list
       will be send only only once per timestep by calling the notify function

       All other topic will be send immediately

       Args:
         port: port subscriber
         host: host subscriber
         notify: notify function
         header: list of all subscribed topis, empty list subscribes all
         num_tasks: number of all tasks
         timestep: timestep used for notify function
    '''
    simple_topic = ['license']
    port_task_step = 3

    def __init__(self, **kwargs):
        '''initialize subscriber
        '''
        threading.Thread.__init__(self)
        self.port = kwargs.get('port', None)
        self.host = kwargs.get('host')
        self.notify = kwargs.get('notify', None)
        self.header = kwargs.get('header')
        self.curve_label = kwargs.get('curve_label', '{}.')
        self.timestep = kwargs.get('timestep', 1)
        self.num_tasks = kwargs.get('num_tasks', 1)
        self.percent_list = [0]*self.num_tasks

        if not self.host:
            self.host = 'localhost'
        if not self.header:
            self.header = [b'']
        self.running = True

        # globals
        self.ylabel_index = 1
        self.notify_send_loop = True if self.notify else False
        self.notify_send_header = set()
        self.notify_send_data = dict()

        # timer function
        if self.notify:
            notify_timerfunc = threading.Timer(0.1, self.send_notify)
            notify_timerfunc.start()

        if b'xyplot' in self.header:
            self.ylabel = self.curve_label.format(self.ylabel_index)
            self.ylabel_index += 1

        context = zmq.Context.instance()
        self.subscriber = context.socket(zmq.SUB)
        from femagtools.multiproc import Engine
        portList = Engine.portPool.portList if Engine.portPool else [self.port]
        if self.num_tasks < len(portList):
            portList = portList[:self.num_tasks]
        for port in portList:
            self.subscriber.connect(f'tcp://{self.host}:{port}')
        logger.debug(f'connect to {portList}')
        self.subscriber.setsockopt(zmq.SUBSCRIBE, self.header[0] if len(self.header) == 1 else b'')
        self.controller = zmq.Context.instance().socket(zmq.PULL)
        self.controller_url = f'inproc://publisher{self.port}'
        try:
            self.controller.bind(self.controller_url)
        except zmq.error.ZMQError:
            pass # ignore

        self.poller = zmq.Poller()
        self.poller.register(self.subscriber, zmq.POLLIN)
        self.poller.register(self.controller, zmq.POLLIN)

    def stop(self):
        ''' stop subscriber
        '''
        socket = zmq.Context.instance().socket(zmq.PUSH)
        socket.connect(self.controller_url)
        socket.send(b"quit")
        socket.close()
        self.running = False

    def send_notify(self):
        ''' timer function to call notify function
        '''
        logger.debug(f"Send loop: {self.notify_send_loop}, timestep: {self.timestep}")
        while self.notify_send_loop:
            if 'progress_logger' in self.notify_send_header:
                # collect data from different threads
                self.notify_send_header.remove('progress_logger')
                numTot = len(self.percent_list)
                d = json.loads(self.notify_send_data.get('progress_logger')[1])
                d['percent'] = sum(self.percent_list) / numTot
                d['subtitle'] = f"{self.percent_list.count(100)} of {numTot}" if numTot > 1 else ''
                self.notify(['progress_logger', json.dumps(d)])
                #logger.debug(f'Send Percent: {d["percent"]} TIMESTEP: {self.timestep}')
            if 'xyplot' in self.notify_send_header:
                self.notify([s.decode('latin1')
                                       for s in self.notify_send_data.get('xyplot')])
                self.notify_send_header.remove('xyplot')

            # simple
            for sdata in SubscriberTask.simple_topic:
                if sdata in self.notify_send_header:
                    self.notify([s.decode('latin1')
                                           for s in self.notify_send_data.get(sdata)])
                    self.notify_send_header.remove(sdata)

            time.sleep(abs(self.timestep))
        logger.debug(f"Send Finished loop: {self.notify_send_loop}")

    def run(self):
        ''' subscriber event loop
        '''
        logger.debug("subscriber is ready, port: %s", {self.port})
        while self.running:
            socks = dict(self.poller.poll())
            if socks.get(self.subscriber) == zmq.POLLIN:
                try:
                    response = self.subscriber.recv_multipart()
                    # Sometimes femag send messages with only len = 1. These messages must be ignored
                    if len(response) < 2:
                        continue
                    # header progress
                    if response[0] == b'progress' and b'progress' in self.header:
                        self.notify_send_header.add('progress_logger')
                        response[0] = b'progress_logger'
                        self.notify_send_data['progress_logger'] = response
                        percent = json.loads(response[1].decode()).get('percent')
                        port = json.loads(response[1].decode()).get('port')
                        hostname = json.loads(response[1].decode()).get('hostname')
                        logger.debug(f"percent: {percent} port: {port}, hostname: {hostname}")
                        self.setMultiProcPercent(percent, hostname, port)
                        continue

                    # header xyplot (add ylabel)
                    if response[0] == b'xyplot' and b'xyplot' in self.header    :
                        d = json.loads(response[1].decode(), strict=False)
                        d['ylabel'] = f"{d.get('ylabel')}_{self.ylabel}" \
                            if d.get('ylabel') else self.ylabel
                        response[1] = json.dumps(d).encode()

                        # timestep negative, immediately update
                        if self.timestep < 0:
                            self.notify([s.decode('latin1') for s in response])
                        else:
                            self.notify_send_data['xyplot'] = response
                            self.notify_send_header.add('xyplot')
                        continue

                    # simple
                    for sdata in SubscriberTask.simple_topic:
                        if response[0] == sdata.encode() and sdata.encode() in self.header:
                            self.notify_send_header.add(sdata)
                            self.notify_send_data[sdata] = response
                            continue

                    if response[0] not in self.header:
                        self.notify([s.decode('latin1') for s in response])

                except Exception:
                    logger.error(
                        "error in subscription message processing", exc_info=True)

            if socks.get(self.controller) == zmq.POLLIN:
                req = self.controller.recv()
                logger.debug("subscriber %s", req)
                break
        self.subscriber.close()
        self.controller.close()
        logger.debug("subscriber stopped")

    def setMultiProcPercent(self, percent, hostname, port):
        ''' set multiproc percent list
        '''
        if not percent:
            return
        from femagtools.multiproc import Engine
        port_list = Engine.portPool.portList if Engine.portPool else [self.port]

        if not port: # old femag
            idx = 0
        elif port_list: # femag publish this port
            idx = port_list.index(port)
        else: # femag publish this port
            logger.info('unknown multiproc progress logic')
            return
        lenPortList = len(port_list)
        lenPercentList = len(self.percent_list)
        if lenPercentList == 1:
            self.percent_list[0] = percent
        else:
            while idx < lenPercentList:
                if self.percent_list[idx] < percent:
                    self.percent_list[idx] = percent
                    break
                elif (idx + lenPortList) < lenPercentList and\
                     self.percent_list[idx + lenPortList] == 0\
                     and percent == 100:
                    break
                else:
                    idx += lenPortList
        #logger.debug(f'{port} percent {percent} index {idx} :: {self.percent_list} len: {lenPortList}')
