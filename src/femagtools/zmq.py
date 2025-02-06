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
    ylabel_index = 1
    curve_label = '{}.'
    # used by static notify func
    percent_list = []
    notify_timerfunc = None
    notify_send_loop = True
    notify = None
    notify_send_header = set()
    notify_send_data = dict()
    # progress, xydata and this entries of this list
    # will send only only once per timestep
    simple_data = ['license']

    def __init__(self, **kwargs):
        threading.Thread.__init__(self)
        context = zmq.Context.instance()
        self.subscriber = context.socket(zmq.SUB)
        self.port = kwargs.get('port', None)
        self.host = kwargs.get('host')
        self.notify = kwargs.get('notify', None)
        SubscriberTask.notify = kwargs.get('notify', None)
        self.header = kwargs.get('header')
        self.num_cur_steps = kwargs.get('num_cur_steps', None)
        SubscriberTask.curve_label = kwargs.get('curve_label', '')
        SubscriberTask.timestep = kwargs.get('timestep', 2)

        if not self.host:
            self.host = 'localhost'
        if not self.header:
            self.header = [b'']
        self.running = True

        # timer function
        if not SubscriberTask.notify_timerfunc:
            SubscriberTask.notify_timerfunc = threading.Timer(0.1, SubscriberTask.send_notify)
            SubscriberTask.notify_send_loop = True
            SubscriberTask.notify_timerfunc.start()

        if b'xyplot' in self.header:
            self.ylabel = self.curve_label.format(SubscriberTask.ylabel_index)
            SubscriberTask.ylabel_index += 1
        if b'progress' in self.header:
            self.protfile = ProtFile(None, self.num_cur_steps)
            self.protId = len(SubscriberTask.percent_list)
            SubscriberTask.percent_list.append(0)  # 0%

        self.subscriber.connect(f'tcp://{self.host}:{self.port}')
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
        self.logger = logger

    def stop(self):
        socket = zmq.Context.instance().socket(zmq.PUSH)
        socket.connect(self.controller_url)
        socket.send(b"quit")
        socket.close()
        self.running = False

    def clear():
        SubscriberTask.ylabel_index = 1
        SubscriberTask.curve_label = '{}.'
        SubscriberTask.notify_timerfunc = None
        SubscriberTask.notify_send_loop = False
        SubscriberTask.notify = None
        SubscriberTask.notify_send_header = set()
        SubscriberTask.notify_send_data = dict()
        SubscriberTask.percent_list = []

    def send_notify():
        logger.debug(f"Send loop: {SubscriberTask.notify_send_loop}, timestep: {SubscriberTask.timestep}")
        while SubscriberTask.notify_send_loop:
            if 'progress_logger' in SubscriberTask.notify_send_header:
                # collect data from different threads
                SubscriberTask.notify_send_header.remove('progress_logger')
                numTot = len(SubscriberTask.percent_list)
                d = json.loads(SubscriberTask.notify_send_data.get('progress_logger')[1])
                d['percent'] = sum(SubscriberTask.percent_list) / numTot
                d['subtitle'] = f"{SubscriberTask.percent_list.count(100)} of {numTot}" if numTot > 1 else ''
                SubscriberTask.notify(['progress_logger', json.dumps(d)])
            if 'xyplot' in SubscriberTask.notify_send_header:
                SubscriberTask.notify([s.decode('latin1')
                                       for s in SubscriberTask.notify_send_data.get('xyplot')])
                SubscriberTask.notify_send_header.remove('xyplot')

            # simple
            for sdata in SubscriberTask.simple_data:
                if sdata in SubscriberTask.notify_send_header:
                    SubscriberTask.notify([s.decode('latin1')
                                           for s in SubscriberTask.notify_send_data.get(sdata)])
                    SubscriberTask.notify_send_header.remove(sdata)

            time.sleep(abs(SubscriberTask.timestep))
        logger.debug(f"Send Finished loop: {SubscriberTask.notify_send_loop}")

    def run(self):
        self.logger.debug("subscriber is ready, port: %s", {self.port})
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
                        SubscriberTask.notify_send_header.add('progress_logger')
                        response[0] = b'progress_logger'
                        SubscriberTask.notify_send_data['progress_logger'] = response
                        SubscriberTask.percent_list[self.protId] = json.loads(response[1].decode()).get('percent')
                        continue

                    # header xyplot (add ylabel)
                    if response[0] == b'xyplot' and b'xyplot' in self.header    :
                        d = json.loads(response[1].decode(), strict=False)
                        d['ylabel'] = f"{d.get('ylabel')}_{self.ylabel}" \
                            if d.get('ylabel') else self.ylabel
                        response[1] = json.dumps(d).encode()

                        # timestep negative, immediately update
                        if SubscriberTask.timestep < 0:
                            self.notify([s.decode('latin1') for s in response])
                        else:
                            SubscriberTask.notify_send_data['xyplot'] = response
                            SubscriberTask.notify_send_header.add('xyplot')
                        continue

                    # simple
                    for sdata in SubscriberTask.simple_data:
                        if response[0] == sdata.encode() and sdata.encode() in self.header:
                            SubscriberTask.notify_send_header.add(sdata)
                            SubscriberTask.notify_send_data[sdata] = response
                            continue

                    if response[0] not in self.header:
                        self.notify([s.decode('latin1') for s in response])

                except Exception:
                    self.logger.error(
                        "error in subscription message processing", exc_info=True)

            if socks.get(self.controller) == zmq.POLLIN:
                req = self.controller.recv()
                self.logger.info("subscriber %s", req)
                break
        self.subscriber.close()
        self.controller.close()
        self.logger.debug("subscriber stopped")
