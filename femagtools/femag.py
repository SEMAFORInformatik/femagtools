# -*- coding: utf-8 -*-
"""
    femagtools.femag
    ~~~~~~~~~~~~~~~~

    Running FEMAG


"""
import subprocess
import os
import glob
import logging
try:
    import zmq
except ImportError:
    pass
import sys
import json
import io
import femagtools.mcv
import femagtools.airgap as ag
import femagtools.fsl
import femagtools.ntib as ntib
import femagtools.config as cfg
import time
import platform
import re
from threading import Thread

logger = logging.getLogger(__name__)


BCHEXT = '.BATCH' if sys.platform.startswith('linux') else '.BCH'  # win32


class FemagError(Exception):
    pass


def subscribe_dev_null(message):
    logger.debug("DevNull: %s", message)
    pass


class BaseFemag(object):
    def __init__(self, workdir, cmd, magnetizingCurves, magnets):
        self.workdir = workdir
        if cmd:
            self.cmd = cmd
        else:
            self.cmd = cfg.get_femag()

        if magnetizingCurves:
            if isinstance(magnetizingCurves,
                          femagtools.mcv.MagnetizingCurve):
                self.magnetizingCurves = magnetizingCurves
            else:
                self.magnetizingCurves = femagtools.mcv.MagnetizingCurve(
                    magnetizingCurves)
        else:
            self.magnetizingCurves = []
        if magnets:
            if isinstance(magnets, femagtools.Magnet):
                self.magnets = magnets
            else:
                self.magnets = femagtools.Magnet(magnets)
        else:
            self.magnets = []

    def copy_magnetizing_curves(self, model, dir=None):
        """extract mc names from model and write files into workdir or dir if given

        Return:
            list of extracted mc names (:obj:`list`)
        """
        dest = dir if dir else self.workdir
        return [self.magnetizingCurves.writefile(m[0], dest, fillfac=m[1])
                for m in model.set_magcurves(
                        self.magnetizingCurves, self.magnets)]

    def create_fsl(self, pmMachine, simulation):
        """create list of fsl commands"""
        model = femagtools.MachineModel(pmMachine)
        self.modelname = model.name
        self.copy_magnetizing_curves(model)

        builder = femagtools.fsl.Builder()
        if simulation:
            return builder.create(model, simulation, self.magnets)
        return builder.create_model(model, self.magnets) + ['save_model("cont")']
        

    def get_log_value(self, pattern, modelname='FEMAG-FSL.log'):
        result = []
        pat = re.compile(pattern)
        with open(os.path.join(self.workdir, 'FEMAG-FSL.log')) as f:
            for line in f:
                m = pat.search(line)
                if m:
                    result.append(float(line.split(':')[-1].split()[0]))
        return result

    def read_bch(self, modelname=None):
        "read most recent BCH/BATCH file and return result"
        # read latest bch file if any
        if not modelname:
            modelname = self._get_modelname_from_log()

        result = femagtools.bch.Reader()
        bchfile_list = sorted(glob.glob(os.path.join(
            self.workdir, modelname+'_[0-9][0-9][0-9]'+BCHEXT)))
        if len(bchfile_list) > 0:
            logger.info("Read BCH {}".format(bchfile_list[-1]))
            with io.open(bchfile_list[-1], encoding='latin1',
                         errors='ignore') as f:
                result.read(f)
        return result

    def read_isa(self, modelname=None):
        "read most recent I7/ISA7 file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()
        import femagtools.isa7
        return femagtools.isa7.read(os.path.join(self.workdir, modelname))

    def read_los(self, modelname=None):
        "read most recent LOS file and return result"
        # read latest los file if any
        if not modelname:
            modelname = self._get_modelname_from_log()

        losfile_list = sorted(glob.glob(os.path.join(
            self.workdir, modelname+'_[0-9][0-9][0-9].LOS')))
        if len(losfile_list) > 0:
            return ntib.read_los(losfile_list[-1])

        return dict()

    def read_airgap_induction(self, modelname=''):
        """read airgap induction"""
        # we need to figure out the number of poles in model
        bch = self.read_bch(modelname)
        return ag.read(os.path.join(self.workdir, 'bag.dat'),
                       bch.machine['p_sim'])

    def _get_modelname_from_log(self):
        """
        Read the modelname from the Femag Log file
        """
        with open(os.path.join(self.workdir, 'FEMAG-FSL.log')) as f:
            for l in f:
                if l.startswith('New model') or l.startswith('Load model'):
                    model = l.split('"')[1]
                    break
        return model


class Femag(BaseFemag):
    """Invoke and control execution of FEMAG

    Args:
        workdir: name of working directory
        cmd: name of femag program (default wfemag64 on windows, xfemag64 on linux)
        magnetizingCurves: collection of lamination material curves
        magnets: collection of magnet material
    """
    def __init__(self, workdir, cmd=None,
                 magnetizingCurves=None, magnets=None):
        super(self.__class__, self).__init__(workdir, cmd,
                                             magnetizingCurves, magnets)

    def run(self, filename, options=['-b'], fsl_args=[]):
        """invoke FEMAG in current workdir

        Args:
            filename: name of file to execute
            options: list of FEMAG options
            fsl_args: list of FSL argument options

        Raises:
            FemagError
        """
        if self.cmd.startswith('wfemag') and \
           '-b' in options and \
           '-m' not in options:
            options.insert(0, '-m')
        args = [self.cmd] + options + [filename] + fsl_args

        basename, ext = os.path.splitext(os.path.basename(filename))
        outname = os.path.join(self.workdir, basename+'.out')
        errname = os.path.join(self.workdir, basename+'.err')
        with open(outname, 'w') as out, open(errname, 'w') as err:
            logger.info('invoking %s', ' '.join(args))
            proc = subprocess.Popen(
                args,
                stdout=out, stderr=err, cwd=self.workdir)

        proc.wait()
        errs = []
        # print femag output
        with io.open(outname, encoding='latin1', errors='ignore') as outfile:
            for l in outfile:
                if l.find('ERROR') > -1:
                    errs.append(l.strip())

        rc = proc.returncode
        logger.info("%s exited with returncode %d (num errs=%d)",
                    self.cmd, rc, len(errs))
        if rc != 0 or errs:
            with io.open(errname, encoding='latin1',
                         errors='ignore') as errfile:
                for l in errfile:
                    errs.append(l.strip())
            errs.insert(0, 'Exit code {}'.format(rc))
            raise FemagError(errs)
        
    def cleanup(self):
        "removes all created files in workdir"
        if not os.path.exists(self.workdir):
            return
        cleanfiles = ('*.B*CH', '*.I*7-*', '*.A*7-*', '*.nc-*',
                      '*.dat', '*.LOS', '*.svg', '*.png', '*.hxy')
        # '*.TMC','*.TMO', '*.PROT'):
        for p in cleanfiles:
            for f in glob.glob(os.path.join(self.workdir, p)):
                os.remove(f)

    def __call__(self, pmMachine, simulation={},
                 options=['-b'], fsl_args=[]):
        """setup fsl file, run calculation and return 
        BCH or LOS results if any."""
        fslfile = 'femag.fsl'
        with open(os.path.join(self.workdir, fslfile), 'w') as f:
            f.write('\n'.join(self.create_fsl(pmMachine,
                                              simulation)))
        if simulation and simulation['calculationMode'] == "pm_sym_loss":
            with open(os.path.join(self.workdir,
                                   self.modelname+'.ntib'), 'w') as f:
                f.write('\n'.join(ntib.create(
                    simulation['speed'],
                    simulation['current'],
                    simulation['angl_i_up'])))
                # TODO: add r1, m

        self.run(fslfile, options, fsl_args)
        if simulation:
            if simulation['calculationMode'] == "pm_sym_loss":
                return self.read_los(self.modelname)
            return self.read_bch(self.modelname)
        return dict(status='ok', message=self.modelname)


class ZmqFemag(BaseFemag):
    """Invoke and control execution of FEMAG with ZeroMQ

    Args:
        port: port number of req socket
        host: hostname (ip addr)
        workdir: name of working directory
        cmd: name of femag program
    """
    def __init__(self, port, host='localhost', workdir='', cmd=None,
                 magnetizingCurves=None, magnets=None):
        super(self.__class__, self).__init__(workdir, cmd,
                                             magnetizingCurves, magnets)
        self.host = host
        self.port = port
        self.ipaddr = ''

        self.request_socket = self.__req_socket()
        self.subscriber_socket = None
        self.proc = None
        self.reader = None

    def close(self):
        if self.reader:
            self.reader.continue_loop = False
        if self.proc:
            self.quit()
        if self.request_socket:
            self.request_socket.close()
            self.request_socket = None
            self.request_socket = self.__req_socket()
        if self.subscriber_socket:
            self.subscriber_socket.close()
            self.subscriber_socket = None

    def __del__(self):
        self.close()
        logger.debug("Destructor ZmqFemag")

    def __req_socket(self):
        """returns a new request client"""
        context = zmq.Context.instance()
        self.request_socket = context.socket(zmq.REQ)
        self.request_socket.connect('tcp://{0}:{1}'.format(
            self.host, self.port))
        #if not self.ipaddr:
        #    if self.host != 'localhost':
        #        inforesp = self.info()
        #        self.ipaddr = json.loads(inforesp[1])['addr']
        #        logger.info("Connected with %s", self.ipaddr)
        #    else:
        #        self.ipaddr = '127.0.0.1'

        return self.request_socket

    def __sub_socket(self):
        """returns a subscriber client"""
        if self.subscriber_socket:
            return self.subscriber_socket
        context = zmq.Context.instance()
        if not self.ipaddr:
            self.ipaddr = '127.0.0.1'
        self.subscriber_socket = context.socket(zmq.SUB)
        self.subscriber_socket.connect(
            'tcp://{0}:{1}'.format(
                self.ipaddr, self.port+1))
        self.subscriber_socket.setsockopt(zmq.SUBSCRIBE, b'')
        self.subscriber_socket.RCVTIMEO = 900  # in milliseconds
        return self.subscriber_socket

    def __is_process_running(self, procId):
        try:
            import psutil
            return psutil.pid_exists(procId)
        except ModuleNotFoundError:
            pass
        # long version, self made
        try:
            if procId > 0:
                if platform.system() == "Windows":
                    #                   if procId in psutil.get_pid_list():
                    proc = subprocess.Popen(["tasklist"],
                                            stdout=subprocess.PIPE)
                    for l in proc.stdout:
                        ls = l.split()
                        try:
                            if str(procId) == ls[1]:
                                return True
                        except IndexError:
                            continue
                else:
                    if not os.kill(procId, 0):
                        return True
        except OSError as e:
            # No such process
            logger.info("OSError: %s", e)
            return False
        except Exception as e:
            # we cannot check processId
            logger.warn("process check error %s", e)
            return True
        return False

    def __is_running(self):
        """check if FEMAG is running in ZMQ mode

        Return:
            True if FEMAG is running, False otherwise

        """
        if not self.request_socket:
            return False

        # call info() and check result
        try:
            ret = [json.loads(s, strict=False) for s in self.info()]
            return ret[0].get('status') == 'ok'
        except Exception:
            pass
        return False

    def send_request(self, msg, pub_consumer=None, timeout=None):
        try:
            # Start the reader thread to get information
            if pub_consumer:
                self.stopStreamReader()
                self.reader = FemagReadStream(self.__sub_socket(), pub_consumer)
                self.reader.setDaemon(True)
                self.reader.start()

            if timeout:
                self.request_socket.setsockopt(zmq.RCVTIMEO, timeout)
                self.request_socket.setsockopt(zmq.LINGER, 0)
            else:
                self.request_socket.setsockopt(zmq.RCVTIMEO, -1)  # default value
                self.request_socket.setsockopt(zmq.LINGER, 30000)  # default value

            while True:
                try:
                    for m in msg[:-1]:
                        self.request_socket.send_string(m, flags=zmq.SNDMORE)
                    if isinstance(msg[-1], list):
                        self.request_socket.send_string('\n'.join(msg[-1]))
                    else:
                        self.request_socket.send_string(msg[-1])
                    return self.request_socket.recv_multipart()
                except zmq.error.Again:
                    pass
                
        except Exception as e:
            #logger.exception("send_request")
            logger.info("send_request: %s Message %s", str(e), msg)
            if timeout:  # only first call raises zmq.error.Again
                logger.info("Femag is not running")
                return [b'{"status":"error", "message":"Femag is not running"}']
            return [b'{"status":"error", "message":"' + str(e).encode() + b'"}']

    def send_fsl(self, fsl, pub_consumer=None, timeout=None):
        """sends FSL commands in ZMQ mode and blocks until commands are processed

        Args:
            fsl: string (or list) of FSL commands
            timeout: The timeout (in milliseconds) to wait for a response
            pub_consumer: callable object to be invoked with publish message 
                      tuple (topic, content)

        Return:
            status
        """
        header = 'FSL'
        logger.debug("Send fsl with fsl: {}, pub_consumer: {}, header: {}".format(
            fsl, pub_consumer, header))

        try:
            # Start the reader thread to get information about the next calculation
            if pub_consumer:
                self.stopStreamReader()
                self.reader = FemagReadStream(self.__sub_socket(), pub_consumer)
                self.reader.setDaemon(True)
                self.reader.start()

            self.request_socket.send_string(header, flags=zmq.SNDMORE)
            logger.debug("Sent header")
            if isinstance(fsl, list):
                self.request_socket.send_string('\n'.join(fsl))
            else:
                self.request_socket.send_string(fsl)
            logger.debug("Sent fsl wait for response")

            if timeout:
                self.request_socket.setsockopt(zmq.RCVTIMEO, timeout)
                self.request_socket.setsockopt(zmq.LINGER, 0)
            else:
                self.request_socket.setsockopt(zmq.RCVTIMEO, -1)  # default value
                self.request_socket.setsockopt(zmq.LINGER, 30000)  # default value
            import datetime
            startTime = datetime.datetime.now() if timeout else None
            while True:
                try:
                    response = self.request_socket.recv_multipart()
                    time.sleep(.5)  # Be sure all messages are arrived over zmq
                    if pub_consumer:
                        self.reader.continue_loop = False
                    return [s.decode('latin1') for s in response]
                except zmq.error.Again as e:
                    logger.info("Again [%s], timeout: %d", str(e), timeout)
                    if startTime:
                        diffTime = datetime.datetime.now() - startTime
                        logger.debug("Diff msec[%d]", diffTime.microseconds)
                        if diffTime.microseconds < timeout:
                            continue
                        else:
                            logger.info("ALERT not running close socket")
                            self.close()
                            return ['{"status":"error", "message":"Femag is not running"}', '{}']

                    continue
        except Exception as e:
            logger.exception("send_fsl")
            logger.error("send_fsl: %s", str(e))
            if timeout:  # only first call raises zmq.error.Again
                return ['{"status":"error", "message":"Femag is not running"}', '{}']
            msg = str(e)
            return ['{"status":"error", "message":"'+msg+'"}', '{}']

    def run(self, options=['-b'], restart=False, procId=None):  # noqa: C901
        """invokes FEMAG in current workdir and returns pid

        Args:
            options: list of FEMAG options
        """
        args = [self.cmd] + options

        if self.__is_running():
            if restart:
                logger.info("must restart")
                self.quit(True)

                # check if process really finished (mq_connection)
                if procId:
                    logger.info("procId: %s", procId)
                    for t in range(200):
                        time.sleep(0.1)
                        if not self.__is_process_running(procId):
                            break
                        logger.info(
                            "femag (pid: '{}') not stopped yet".format(
                                procId))
                    logger.info("Stopped procId: %s", procId)
            else:
                try:
                    with open(os.path.join(self.workdir,
                                           'femag.pid'), 'r') as pidfile:
                        procId = int(pidfile.readline())
                except Exception:
                    pass
                return procId

        basename = str(self.port)
        args.append(basename)

        outname = os.path.join(self.workdir, basename+'.out')
        errname = os.path.join(self.workdir, basename+'.err')
        with open(outname, 'w') as out, open(errname, 'w') as err:
            logger.debug('invoking %s', ' '.join(args))
            self.proc = subprocess.Popen(
                args,
                stdout=out, stderr=err, cwd=self.workdir)
            rc = self.proc.poll()
            if rc:
                raise RuntimeError('Process "{}" exited with {}'.format(
                    ' '.join(args), rc))

        # check if mq is ready for listening
        lcount = 10
        for t in range(lcount):
            time.sleep(0.1)
            if self.__is_running():
                logger.info("femag (pid: '{}') is listening".format(
                    self.proc.pid))
                break
            else:
                # reopen request socket
                logger.debug('reopen request port')
                self.request_socket.close()
                self.request_socket = None
                if t == (lcount-1):
                    # 10 attempts fails, abort
                    logger.info('abort (starting femag)')
                    self.proc.wait(1000)
                    rc = self.proc.returncode
                    self.proc = None
                    return rc  # self.proc.returncode
                self.request_socket = self.__req_socket()

        # write femag.pid
        with open(os.path.join(self.workdir, 'femag.pid'), 'w') as pidfile:
            pidfile.write("{}\n".format(self.proc.pid))
        return self.proc.pid

    def quit(self, save_model=False):
        """terminates femag"""

        if not self.__is_running():
            logger.info("Femag already stopped")
            if self.proc:
                self.proc.wait()
                self.proc = None
            return

        # send exit flags
        f = '\n'.join(['exit_on_end = true',
                       'exit_on_error = true'])
        response = self.send_fsl(f, pub_consumer=subscribe_dev_null)

        # send quit command
        try:
            response = [r.decode('latin1')
                        for r in self.send_request(
                                ['CONTROL', 'quit'], timeout=10000,
                                pub_consumer=subscribe_dev_null)]
        except Exception as e:
            logger.error("Femag Quit zmq message %s", e)

        logger.debug("Sent QUIT to femag %s", response)
        # if query, send a answer
        obj = json.loads(response[0])
        logger.debug("status: {}".format(obj['status']))
        if obj['status'] == 'Query':
            logger.info('query: %s => %s',
                        obj['message'], 'saved' if save_model else 'not saved')

            # Only send one msg
            response = self.request_socket.send_string(
                'Ok' if save_model else 'Cancel')
        
        if self.proc:
            self.proc.wait()
            self.proc = None
        return response

    def upload(self, files):
        """upload file or files
        returns number of transferred bytes
        (FEMAG 8.5 Rev 3282 or greater only)
        """
        ret = []
        total = 0
        chunk_size = 20*1024
        fnames = files
        if isinstance(files, str):
            fnames = [files]
        for fn in fnames:
            self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
            basename = os.path.basename(fn)
            logger.info("upload %s --> %s", fn, basename)
            self.request_socket.send_string('upload = {}'.format(basename),
                                            flags=zmq.SNDMORE)
            with open(fn, mode="rb") as file:
                while True:
                    data = file.read(chunk_size)
                    if not data:
                        break
                    more = 0 if len(data) < chunk_size else zmq.SNDMORE
                    self.request_socket.send(data, flags=more)
                    total += len(data)
            ret += [s.decode('latin1')
                    for s in self.request_socket.recv_multipart()]
        return ret
    
    def cleanup(self):
        """remove all FEMAG files in working directory 
        (FEMAG 8.5 Rev 3282 or greater only)"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'cleanup'])]
    
    def release(self):
        """signal finish calculation task to load balancer to free resources
        (Docker Cloud environment only)
        """
        return [r.decode('latin1')
                for r in self.send_request(['close'])]

    def info(self, timeout=2000):
        """get various resource information 
        (FEMAG 8.5 Rev 3282 or greater only)"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'info'],
                                           timeout=timeout)]

    def publishLevel(self, level):
        """set publish level"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'publish = {}'
                                            .format(level)], timeout=10000)]

    def getfile(self, filename=''):
        """get file (FEMAG 8.5 Rev 3282 or greater only)"""
        response = self.send_request(
            ['CONTROL', 'getfile = {}'.format(filename)])
        return [response[0].decode('latin1'),
                response[1] if len(response) else b'']
        
    def exportsvg(self, fslcmds):
        """get svg format from fsl commands (if any graphic created)
        (since FEMAG 8.5 Rev 3343) """
        response = self.send_request(['SVG', fslcmds])
        rc = json.loads(response[0].decode('latin1'))
        if rc['status'] == 'ok':
            return self.getfile(rc['result_file'][0])
        return [s.decode('latin1') for s in response]
            
    def stopStreamReader(self):
        if self.reader:
            logger.debug("stop stream reader")
            self.reader.continue_loop = False
            self.reader.join()

    def copy_magnetizing_curves(self, model, dir=None):
        """extract mc names from model and write files into workdir or dir if given
           and upload to Femag

        Return:
            list of extracted mc names (:obj:`list`)
        """
        dest = dir if dir else self.workdir
        for m in model.set_magcurves(
                self.magnetizingCurves, self.magnets):
            f = self.magnetizingCurves.writefile(m[0], dest, fillfac=m[1])
            self.upload(os.path.join(dest, f))
    
    def __call__(self, pmMachine, simulation, pub_consumer=None):
        """setup fsl file, run calculation and return BCH results
        Args:
          pmMachine: dict with machine parameters or name of model
          simulation; dict with simulation parameters

        Raises:
           FemagError
        """
        if isinstance(pmMachine, str):
            modelpars = dict(name=pmMachine)
        else:
            modelpars = pmMachine
        if 'exit_on_end' not in modelpars:
            modelpars['exit_on_end'] = 'false'
        if 'exit_on_error' not in modelpars:
            modelpars['exit_on_error'] = 'false'
        response = self.send_fsl(['save_model("close")'] +
                                 self.create_fsl(modelpars,
                                                 simulation),
                                 pub_consumer=pub_consumer)
        r = json.loads(response[0])
        if r['status'] != 'ok':
            raise FemagError(r['message'])

        result_file = r['result_file'][0]
        if simulation['calculationMode'] == "pm_sym_loss":
            return self.read_los(self.modelname)

        status, content = self.getfile(result_file)
        r = json.loads(status)
        if r['status'] == 'ok':
            bch = femagtools.bch.Reader()
            return bch.read(content.decode('latin1'))
        raise FemagError(r['message'])


class FemagReadStream(Thread):
    def __init__(self, sub_socket, pub_consumer):
        Thread.__init__(self)
        logger.debug("Initialize reader thread")
        self.sub_socket = sub_socket
        # TODO: use a function to handle calls without pub_consumer
        self.pub_consumer = pub_consumer if pub_consumer else self.dummy_pub_consumer
        self.continue_loop = True

    def dummy_pub_consumer(self, data):
        """This dummy method is used when no pub_consumer is defined.
        """
        logger.warn(
            "No pub_consumer defined, fallback to dummy\n >{}".format(data))

    def run(self):
        """Listen for messages from femag as long as the continue_loop is True
        the sub_socket has a timeout to finish the loop and thread after
        the continue_loop is set to False and no messages are arrived.
        """
        logger.debug(
            "Start thread with while condition: {}".format(self.continue_loop))
        while self.continue_loop:
            try:
                response = self.sub_socket.recv_multipart()
                # Sometimes femag send messages with only len = 1. These messages must be ignored
                if len(response) < 2:
                    continue
                # Call the pub_consumer function
                self.pub_consumer([s.decode('latin1')
                                   for s in response])
            # The subscriber_socket has a timeout of 900 mil sec. If no answer was arrived
            # this exception is raised - Ignore
            except zmq.error.Again:
                continue
            # Any other exception is shown in the error log
            except Exception as e:
                logger.error("error in reading output from femag: {}".format(e))
                continue
        logger.debug("Exit reader thread")
