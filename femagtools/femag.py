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
        return [self.magnetizingCurves.writefile(m, dest)
                for m in model.set_magcurves(
                        self.magnetizingCurves, self.magnets)]

    def create_fsl(self, pmMachine, operatingConditions):
        """create list of fsl commands"""
        model = femagtools.MachineModel(pmMachine)
        self.modelname = model.name
        self.copy_magnetizing_curves(model)

        builder = femagtools.fsl.Builder()
        return builder.create(model, operatingConditions, self.magnets)

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

    def read_airgap_induction(self):
        """read airgap induction"""
        return ag.read(os.path.join(self.workdir, 'bag.dat'))
    
    def _get_modelname_from_log(self):
        """
        Read the modelname from the Femag Log file
        """
        with open(os.path.join(self.workdir, 'FEMAG-FSL.log')) as f:
            for l in f:
                if l.startswith('New model') or l.startswith('Load model'):
                    model = l.split()[2].replace('"', '').replace(',', '')
                    break
        return model


class Femag(BaseFemag):
    """Invoke and control execution of FEMAG

    Args:
        workdir: name of working directory
        cmd: name of femag program
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
        print(outname)
        with io.open(outname, encoding='latin1', errors='ignore') as outfile:
            for l in outfile:
                print(l.strip())
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
        cleanfiles = ('*.B*CH', '*.I*7-*', '*.A*7-*',
                      '*.dat', '*.LOS', '*.svg', '*.png')
        # '*.TMC','*.TMO', '*.PROT', '*.hxy'):
        for p in cleanfiles:
            for f in glob.glob(os.path.join(self.workdir, p)):
                os.remove(f)

    def __call__(self, pmMachine, operatingConditions,
                 options=['-b'], fsl_args=[]):
        """setup fsl file, run calculation and return BCH results"""
        fslfile = 'femag.fsl'
        with open(os.path.join(self.workdir, fslfile), 'w') as f:
            f.write('\n'.join(self.create_fsl(pmMachine,
                                              operatingConditions)))
        if operatingConditions['calculationMode'] == "pm_sym_loss":
            with open(os.path.join(self.workdir,
                                   self.modelname+'.ntib'), 'w') as f:
                f.write('\n'.join(ntib.create(
                    operatingConditions['speed'],
                    operatingConditions['current'],
                    operatingConditions['angl_i_up'])))
                # TODO: add r1, m

        self.run(fslfile, options, fsl_args)
        if operatingConditions['calculationMode'] == "pm_sym_loss":
            return self.read_los(self.modelname)
        return self.read_bch(self.modelname)


class ZmqFemag(BaseFemag):
    """Invoke and control execution of FEMAG with ZeroMQ

    Args:
        workdir: name of working directory
        cmd: name of femag program
    """
    def __init__(self, workdir, port, host='localhost', cmd=None,
                 magnetizingCurves=None, magnets=None):
        super(self.__class__, self).__init__(workdir, cmd,
                                             magnetizingCurves, magnets)
        self.host = host
        self.port = port
        self.request_socket = None
        self.subscriber_socket = None
        self.proc = None
        self.reader = None

    def __del__(self):
        if self.reader:
            self.reader.continue_loop = False
        logger.debug("Destructor ZmqFemag")

    def __req_socket(self):
        """returns a new request client"""
        if self.request_socket:
            return self.request_socket
        context = zmq.Context.instance()
        self.request_socket = context.socket(zmq.REQ)
        self.request_socket.connect('tcp://{0}:{1}'.format(
            self.host, self.port))
        self.request_socket.setsockopt(zmq.LINGER, 500)
        return self.request_socket

    def __sub_socket(self):
        """returns a subscriber client"""
        if self.subscriber_socket:
            return self.subscriber_socket
        context = zmq.Context.instance()
        self.subscriber_socket = context.socket(zmq.SUB)
        self.subscriber_socket.connect(
            'tcp://{0}:{1}'.format(
                self.host, self.port+1))
        self.subscriber_socket.setsockopt(zmq.SUBSCRIBE, b'')
        self.subscriber_socket.RCVTIMEO = 900  # in milliseconds
        return self.subscriber_socket

    def __is_process_running(self, procId):
        try:
            import psutil
            return psutil.pid_exists(procId)
        except:
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
                        except:
                            continue
                else:
                    if not os.kill(procId, 0):
                        return True
        except OSError as e:
            # No such process
            logger.info("OSError: '{}'\n".format(str(e)))
            return False
        except Exception:
            # we cannot check processId
            logger.info("Error: unknown\n")
            return True
        return False

    def __is_running(self, timeout=1500):
        """check if FEMAG is running in ZMQ mode

        Args:
            timeout: The timeout (in milliseconds) to wait for a response

        Return:
            True if FEMAG is running, False otherwise

        """
        try:
            request_socket = self.__req_socket()
            request_socket.send_string("FSL", flags=zmq.SNDMORE)
            request_socket.send_string("testvar=0")
            poller = zmq.Poller()
            # use POLLIN for recv, POLLOUT for send
            poller.register(request_socket, zmq.POLLIN)
            if poller.poll(timeout):  # ok, femag is running
                logger.info('femag is running for %s',
                            self.request_socket.getsockopt_string(
                                zmq.IDENTITY))
                request_socket.recv_multipart()
                return True
        except Exception as e:
            logger.error(e)
        self.request_socket.close()
        self.request_socket = None
        logger.debug("femag is not running")
        return False

    def send_fsl(self, fsl, callback=None, header='FSL', timeout=None):
        """sends FSL commands in ZMQ mode and blocks until commands are processed

        Args:
            fsl: string of FSL commands
            timeout: The timeout (in milliseconds) to wait for a response

        Return:
            response
        """
        logger.debug("Send fsl with fsl: {}, callback: {}, header: {}".format(
            fsl, callback, header))

        try:
            # Start the reader thread to get information about the next calculation
            if callback:
                self.stopStreamReader()
                self.reader = FemagReadStream(self.__sub_socket(), callback)
                self.reader.setDaemon(True)
                self.reader.start()

            request_socket = self.__req_socket()
            request_socket.send_string(header, flags=zmq.SNDMORE)
            logger.debug("Sent header")
            request_socket.send_string(fsl)
            logger.debug("Sent fsl wait for response")

            if timeout:
                request_socket.setsockopt(zmq.RCVTIMEO, timeout)
                request_socket.setsockopt(zmq.LINGER, 0)
            else:
                request_socket.setsockopt(zmq.RCVTIMEO, -1)  # default value
                request_socket.setsockopt(zmq.LINGER, 30000)  # default value
            try:
                response = request_socket.recv_multipart()
                time.sleep(.5)  # Be sure all messages are arrived over zmq
                if callback:
                    self.reader.continue_loop = False
                return [s.decode() for s in response]
            except zmq.error.Again as e:
                logger.info("Again [%s], timeout: %d", str(e), timeout)
                return ['{"status":"error", "message":"Femag is not running"}', '{}']
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

        Raises:
            FemagError
        """
        args = [self.cmd] + options

        if self.__is_running():
            if restart:
                logging.info("must restart")
                self.quit(True)

                # check if process really finished (mq_connection)
                logger.info("procId: %s", procId)
                if procId:
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

        if self.request_socket:
            self.request_socket.close()
        self.request_socket = None

        basename = str(self.port)
        args.append(basename)

        outname = os.path.join(self.workdir, basename+'.out')
        errname = os.path.join(self.workdir, basename+'.err')
        with open(outname, 'w') as out, open(errname, 'w') as err:
            logger.info('invoking %s', ' '.join(args))
            self.proc = subprocess.Popen(
                args,
                stdout=out, stderr=err, cwd=self.workdir)

        # check if mq is ready for listening
        for t in range(200):
            time.sleep(0.1)
            if self.__is_running():
                logger.info("femag (pid: '{}') is listening".format(
                    self.proc.pid))
                break

        # write femag.pid
        logger.info("ready %s", self.proc.pid)
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
        self.__req_socket()
        f = '\n'.join(['exit_on_end = true',
                       'exit_on_error = true'])
        response = self.send_fsl(f)

        # send quit command
        try:
            response = self.send_fsl('quit', header='CONTROL')
#                                     recvflags=zmq.NOBLOCK)
        except Exception as e:
            logger.error("Femag Quit zmq message %s", e)

        logger.info("Sent QUIT to femag %s", response)
        # if query, send a answer
        obj = json.loads(response[0])
        logger.debug("status: {}".format(obj['status']))
        if obj['status'] == 'Query':
            logger.info('query: %s => %s',
                        obj['message'], 'saved' if save_model else 'not saved')

            # Only send one msg
            request_socket = self.__req_socket()
            response = request_socket.send_string(
                'Ok' if save_model else 'Cancel')
        
        if self.proc:
            self.proc.wait()
            self.proc = None
        return response

    def stopStreamReader(self):
        if self.reader:
            logger.debug("stop stream reader")
            self.reader.continue_loop = False
            self.reader.join()

    def __call__(self, pmMachine, operatingConditions):
        """setup fsl file, run calculation and return BCH results"""
        self.send_fsl('\n'.join(self.create_fsl(pmMachine,
                                                operatingConditions)))
        return self.read_bch(self.modelname)


class FemagReadStream(Thread):
    def __init__(self, sub_socket, callback):
        Thread.__init__(self)
        logger.debug("Initialize reader thread")
        self.sub_socket = sub_socket
        # TODO: use a function to handle calls without callback
        self.callback = callback if callback else self.dummy_callback
        self.continue_loop = True

    def dummy_callback(self, data):
        """This dummy method is used when no callback is defined.
        """
        logger.warn(
            "No callback defined, fallback to dummy\n >{}".format(data))

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
                # Call the callback function
                self.callback([s.decode() for s in response])
            # The subscriber_socket has a timeout of 900 mil sec. If no answer was arrived
            # this exception is raised - Ingnore
            except zmq.error.Again:
                continue
            # Any other exception is shown in the error log
            except Exception as e:
                logger.error("error in reading output from femag {}".format(e))
                continue
        logger.debug("Exit reader thread")
