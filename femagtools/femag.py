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
import femagtools.model
import femagtools.magnet
import femagtools.conductor
import femagtools.windings
import femagtools.mcv
import femagtools.airgap as ag
import femagtools.fsl
import femagtools.ntib as ntib
import femagtools.config as cfg
import time
import re
import threading

logger = logging.getLogger(__name__)


BCHEXT = 'BATCH' if sys.platform.startswith('linux') else 'BCH'  # win32


def get_shortCircuit_parameters(bch, nload):
    try:
        if nload < 0:
            nload = 0
        if nload > 2:
            nload = 2
        if nload > 0:
            dqld = bch.dqPar['ld']
            dqlq = bch.dqPar['lq']
            dqpsim = bch.dqPar['psim']
            if len(dqld) <= nload or len(dqlq) <= nload or len(dqpsim) <= nload:
                ld = dqld[-1]/bch.armatureLength
                lq = dqlq[-1]/bch.armatureLength
                psim = dqpsim[-1]/bch.armatureLength
            else:
                ld = dqld[nload-1]/bch.armatureLength
                lq = dqlq[nload-1]/bch.armatureLength
                psim = dqpsim[nload-1]/bch.armatureLength
        else:
            ld = bch.machine['ld']/bch.armatureLength
            lq = bch.machine['lq']/bch.armatureLength
            psim = bch.machine['psim']/bch.armatureLength
        return dict(
            r1=bch.machine['r1'],
            ld=ld,
            lq=lq,
            psim=psim,
            num_pol_pair=bch.machine['p'],
            fc_radius=bch.machine['fc_radius'],
            lfe=bch.armatureLength/1e3,
            pocfilename=bch.machine['pocfile'],
            num_par_wdgs=bch.machine['num_par_wdgs'],
            calculationMode='shortcircuit')
    except (KeyError, AttributeError, IndexError):
        raise FemagError("missing pm/Rel-Sim results")


class FemagError(Exception):
    pass


class BaseFemag(object):
    def __init__(self, workdir, cmd, magnetizingCurves, magnets, condMat):
        self.workdir = workdir
        self.magnets = []
        self.magnetizingCurves = []
        self.condMat = []
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

        if magnets:
            if isinstance(magnets, femagtools.magnet.Magnet):
                self.magnets = magnets
            else:
                self.magnets = femagtools.magnet.Magnet(magnets)

        if condMat:
            if isinstance(condMat, femagtools.conductor.Conductor):
                self.condMat = condMat
            else:
                self.condMat = femagtools.conductor.Conductor(condMat)

    def copy_winding_file(self, name, wdg):
        wdg.write(name, self.workdir)

    def copy_magnetizing_curves(self, model, dir=None):
        """extract mc names from model and write files into workdir or dir if given

        Return:
            list of extracted mc names (:obj:`list`)
        """
        dest = dir if dir else self.workdir
        return [self.magnetizingCurves.writefile(m[0], dest, fillfac=m[1])
                for m in model.set_magcurves(
            self.magnetizingCurves, self.magnets)]

    def create_wdg_def(self, model):
        name = 'winding'
        w = femagtools.windings.Winding(
            dict(
                Q=model.stator['num_slots'],
                p=model.poles//2,
                m=len(model.windings['wdgdef']),
                windings=model.windings['wdgdef']))
        self.copy_winding_file(name, w)
        return name

    def create_fsl(self, machine, simulation):
        """create list of fsl commands"""
        self.model = femagtools.model.MachineModel(machine)
        self.modelname = self.model.name
        self.copy_magnetizing_curves(self.model)
        try:
            if 'wdgdef' in self.model.windings:
                self.model.windings['wdgfile'] = self.create_wdg_def(
                    self.model)
        except AttributeError:
            pass
        builder = femagtools.fsl.Builder()
        if simulation:
            return builder.create(self.model, simulation, self.magnets, self.condMat)
        return builder.create_model(self.model,
                                    self.magnets,
                                    self.condMat) + ['save_model("cont")']

    def get_log_value(self, pattern, modelname='FEMAG-FSL.log'):
        result = []
        pat = re.compile(pattern)
        with open(os.path.join(self.workdir, 'FEMAG-FSL.log')) as f:
            for line in f:
                m = pat.search(line)
                if m:
                    result.append(float(line.split(':')[-1].split()[0]))
        return result

    def get_bch_file(self, modelname, offset=0):
        return self.get_result_file(modelname, BCHEXT, offset)

    def get_asm_file(self, modelname, offset=0):
        return self.get_result_file(modelname, 'ASM', offset)

    def get_result_file(self, modelname, ext, offset=0):
        """return latest result (bch, asm) file (if any)"""
        filelist = sorted(glob.glob(os.path.join(
            self.workdir, modelname+'_[0-9][0-9][0-9].'+ext)))
        if(filelist):
            return filelist[-1-offset]
        return ''

    def read_asm(self, modelname=None, offset=0):
        "read most recent ASM file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()

        asmfile = self.get_asm_file(modelname, offset)
        if asmfile:
            logger.info("Read ASM {}".format(asmfile))
            return femagtools.asm.read(asmfile)
        return {}

    def read_bch(self, modelname=None, offset=0):
        "read most recent BCH/BATCH file and return result"
        # read latest bch file if any
        if not modelname:
            modelname = self._get_modelname_from_log()

        result = femagtools.bch.Reader()
        bchfile = self.get_bch_file(modelname, offset)
        if bchfile:
            logger.info("Read BCH {}".format(bchfile))
            with io.open(bchfile, encoding='latin1',
                         errors='ignore') as f:
                result.read(f)
        return result

    def read_isa(self, modelname=None):
        "read most recent I7/ISA7 file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()
        import femagtools.isa7
        return femagtools.isa7.read(os.path.join(self.workdir, modelname))

    def read_nc(self, modelname=None):
        "read most recent NC file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()
        import femagtools.nc
        return femagtools.nc.read(os.path.join(self.workdir, modelname))

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

    def read_airgap_induction(self, modelname='', offset=0):
        """read airgap induction"""
        # we need to figure out the number of poles in model
        bch = self.read_bch(modelname, offset)
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
        condMat: collection of conductor material
    """

    def __init__(self, workdir, cmd=None,
                 magnetizingCurves=None, magnets=None, condMat=[]):
        super(self.__class__, self).__init__(workdir, cmd,
                                             magnetizingCurves, magnets, condMat)

    def run(self, filename, options=['-b'], fsl_args=[]):
        """invoke FEMAG in current workdir

        Args:
            filename: name of file to execute
            options: list of FEMAG options
            fsl_args: list of FSL argument options

        Raises:
            FemagError
        """
        if self.cmd.find('wfemag') > -1 and \
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

    def __call__(self, machine, simulation={},
                 options=['-b'], fsl_args=[]):
        """setup fsl file, run calculation and return
        BCH or LOS results if any."""
        fslfile = 'femag.fsl'
        with open(os.path.join(self.workdir, fslfile), 'w') as f:
            f.write('\n'.join(self.create_fsl(machine,
                                              simulation)))
        if simulation:
            if 'poc' in simulation:
                with open(os.path.join(self.workdir,
                                       simulation['pocfilename']), 'w') as f:
                    f.write('\n'.join(simulation['poc'].content()))
            if simulation['calculationMode'] == "pm_sym_loss":
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

            if simulation['calculationMode'] == 'asyn_motor':
                return self.read_asm(self.modelname)

            bch = self.read_bch(self.modelname)
            if simulation['calculationMode'] == 'pm_sym_fast':
                if simulation.get('shortCircuit', False):
                    logger.info("short circuit simulation")
                    simulation.update(
                        get_shortCircuit_parameters(bch,
                                                    simulation.get('initial', 2)))
                    builder = femagtools.fsl.Builder()
                    fslcmds = (builder.open_model(self.model) +
                               builder.create_shortcircuit(simulation))
                    with open(os.path.join(self.workdir, fslfile), 'w') as f:
                        f.write('\n'.join(fslcmds))
                    self.run(fslfile, options)
                    bchfile = self.get_bch_file(self.modelname)
                    if bchfile:
                        bchsc = femagtools.bch.Reader()
                        logger.info("Read BCH {}".format(bchfile))
                        with io.open(bchfile, encoding='latin1',
                                     errors='ignore') as f:
                            bchsc.read(f)
                    bch.scData = bchsc.scData
                    for w in bch.flux:
                        try:
                            bch.flux[w] += bchsc.flux[w]
                            bch.flux_fft[w] += bchsc.flux_fft[w]
                        except (KeyError, IndexError):
                            logging.debug(
                                "No additional flux data in sc simulation")
                            break

                    bch.torque += bchsc.torque
                    bch.demag += bchsc.demag
            return bch
        return dict(status='ok', message=self.modelname)


class FemagTask(threading.Thread):
    def __init__(self, port, args, workdir, logdir):
        threading.Thread.__init__(self)
        self.port = port
        self.args = args + [str(self.port)]
        self.returncode = None
        self.workdir = workdir
        self.logdir = logdir

    def run(self):
        logger.info("femag is ready on port %d workdir %s",
                    self.port, self.workdir)
        outname = os.path.join(self.logdir, f'femag-{self.port}.out')
        errname = os.path.join(self.logdir, f'femag-{self.port}.err')
        with open(outname, 'w') as out, open(errname, 'w') as err:
            self.proc = subprocess.Popen(
                self.args,
                stdout=out, stderr=err, cwd=self.workdir)

        self.returncode = self.proc.wait()


class SubscriberTask(threading.Thread):
    def __init__(self, port, host, notify):
        threading.Thread.__init__(self)
        context = zmq.Context.instance()
        self.subscriber = context.socket(zmq.SUB)
        if not host:
            host = 'localhost'
        self.subscriber.connect(f'tcp://{host}:{port}')
        self.subscriber.setsockopt(zmq.SUBSCRIBE, b'')
        self.controller = zmq.Context.instance().socket(zmq.PULL)
        self.controller_url = 'inproc://publisher'
        self.controller.bind(self.controller_url)
        self.poller = zmq.Poller()
        self.poller.register(self.subscriber, zmq.POLLIN)
        self.poller.register(self.controller, zmq.POLLIN)
        self.logger = logger
        self.notify = notify

    def stop(self):
        socket = zmq.Context.instance().socket(zmq.PUSH)
        socket.connect(self.controller_url)
        socket.send(b"quit")
        socket.close()

    def run(self):
        self.logger.info("subscriber is ready")
        while True:
            socks = dict(self.poller.poll())
            if socks.get(self.subscriber) == zmq.POLLIN:
                try:
                    response = self.subscriber.recv_multipart()
                    # Sometimes femag send messages with only len = 1. These messages must be ignored
                    if len(response) < 2:
                        continue
                    self.notify([s.decode('latin1') for s in response])

                except Exception:
                    self.logger.error(
                        "error in subscription message processing", exc_info=True)

            if socks.get(self.controller) == zmq.POLLIN:
                req = self.controller.recv()
                self.logger.info(req)
                break

        self.logger.debug("subscriber stopped")


class ZmqFemag(BaseFemag):
    """Invoke and control execution of FEMAG with ZeroMQ

    Args:
        port: port number of req socket
        host: hostname (ip addr)
        workdir: name of working directory
        logdir: name of logging directory (default is workdir/log)
        cmd: name of femag program
    """

    def __init__(self, port, host='localhost', workdir='', logdir='',
                 cmd=None,
                 magnetizingCurves=None, magnets=None, condMat=[]):
        super(self.__class__, self).__init__(workdir, cmd,
                                             magnetizingCurves, magnets, condMat)
        self.host = host
        self.port = port
        self.femaghost = ''
        self.logdir = logdir if logdir else os.path.join(workdir, 'log')
        if not os.path.exists(self.logdir):
            os.makedirs(self.logdir)
        self.request_socket = None
        self.subscriber = None

    def close(self):
        if self.subscriber:
            self.subscriber.stop()

        if self.request_socket:
            logger.debug("close request_socket")
            self.request_socket.close()
            self.request_socket = None

        logger.debug("done")

    def __del__(self):
        self.close()
        logger.debug("Destructor ZmqFemag")

    def stopFemag(self):
        if self.femagTask and self.femagTask.proc:
            logger.info("stopFemagTask")
            self.femagTask.proc.kill()
            self.femagTask = None
            self.request_socket = None
            self.subscriber.stop()
            self.subscriber = None
            logger.info("stopFemagTask Done")
        else:
            logger.warning("stopFemag not implemented")

    def __req_socket(self):
        """returns a new request client"""
        context = zmq.Context.instance()
        request_socket = context.socket(zmq.REQ)
        url = 'tcp://{0}:{1}'.format(
            self.host, self.port)
        logger.debug("connect req socket %s", url)
        request_socket.connect(url)
        return request_socket

    def subscribe(self, notify):
        """attaches a notify function"""
        logger.info("Subscribe on '%s'", self.femaghost)
        if self.subscriber is None:
            self.subscriber = SubscriberTask(
                self.port+1, self.femaghost, notify)
            self.subscriber.start()
        else:
            # reattach?
            self.subscriber.notify = notify
        return

    def __is_running(self):
        try:
            return self.femagTask.proc and self.femagTask.returncode is None
        except:
            pass
        return False

    def send_request(self, msg, pub_consumer=None, timeout=None):
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        if not msg:
            return [b'{"status":"ignored",{}}']
        if timeout:
            self.request_socket.setsockopt(zmq.RCVTIMEO, timeout)
            self.request_socket.setsockopt(zmq.LINGER, 0)
        else:
            self.request_socket.setsockopt(zmq.RCVTIMEO, -1)  # default value

        for m in msg[:-1]:
            self.request_socket.send_string(m, flags=zmq.SNDMORE)
        if isinstance(msg[-1], list):
            self.request_socket.send_string('\n'.join(msg[-1]))
        else:
            self.request_socket.send_string(msg[-1])

        errmsg = ''
        while True:
            try:
                return self.request_socket.recv_multipart()
            except zmq.error.Again as e:
                # logger.exception("send_request")
                errmsg = str(e)
                logger.warning("send_request: %s Message %s", str(e), msg)
                continue
        logger.info("oops")
        return [b'{"status":"error", "message":"' + errmsg.encode() + b'"}']

    def send_fsl(self, fsl, timeout=None):
        """sends FSL commands in ZMQ mode and blocks until commands are processed

        Args:
            fsl: string (or list) of FSL commands
            timeout: The maximum time (in milliseconds) to wait for a response

        Return:
            status
        """
        header = 'FSL'
        self.request_socket.send_string(header, flags=zmq.SNDMORE)

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
        import datetime
        startTime = datetime.datetime.now() if timeout else None
        while True:
            try:
                response = self.request_socket.recv_multipart()
                # NOTE: femag encoding is old school
                return [s.decode('latin1') for s in response]
            except zmq.error.Again as e:
                logger.info("Again [%s], timeout: %d", str(e), timeout)
                if not startTime:
                    continue

                diffTime = datetime.datetime.now() - startTime
                logger.debug("Diff msec[%d]", diffTime.microseconds)
                if diffTime.microseconds < timeout:
                    continue

                logger.info("ALERT not running close socket")
                self.close()
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
        """
        if self.__is_running():
            if restart:
                logger.info("must restart")
                self.quit(True)
                self.femagTask.join()
                logger.info("Stopped procId: %s", self.femagTask.proc.pid)
            else:
                return self.femagTask.proc.pid

        if self.cmd.find('wfemag') > -1 and \
           '-b' in options and \
           '-m' not in options:
            options.insert(0, '-m')

        args = [self.cmd] + options
        self.femagTask = FemagTask(self.port, args, self.workdir, self.logdir)
        self.femagTask.start()
        if not self.request_socket:
            self.request_socket = self.__req_socket()

        # check if mq is ready for listening
        lcount = 10
        for t in range(lcount):
            time.sleep(0.1)
            if self.__is_running():
                if self.femagTask.proc:
                    logger.info("femag (pid: '{}') is listening".format(
                        self.femagTask.proc.pid))
                break

        return self.femagTask.proc.pid

    def quit(self, save_model=False):
        """terminates femag"""

        if not self.__is_running():
            logger.info("Femag already stopped")
            return

        if self.subscriber:
            self.subscriber.stop()
            self.subscriber = None

        # send exit flags
        f = '\n'.join(['exit_on_end = true',
                       'exit_on_error = true'])
        response = self.send_fsl(f)

        # send quit command
        try:
            response = [r.decode('latin1')
                        for r in self.send_request(
                ['CONTROL', 'quit'], timeout=2000)]
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

    def copyfile(self, filename, dirname):
        """copy filename to dirname (FEMAG 9.2)"""
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
        self.request_socket.send_string(f'copyfile {filename} {dirname}')
        return [r.decode() for r in self.request_socket.recv_multipart()]

    def change_case(self, dirname):
        """change case to dirname (FEMAG 9.2)"""
        logger.info("change_case to :  %s", dirname)
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
        self.request_socket.send_string(f'casedir {dirname}')
        return [r.decode() for r in self.request_socket.recv_multipart()]

    def delete_case(self, dirname):
        """delete case dir (FEMAG 9.2)"""
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
        self.request_socket.send_string(f'casedir -{dirname}')
        return [r.decode() for r in self.request_socket.recv_multipart()]

    def clear(self, timeout=2000):
        """clear lua script session"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'clear'],
                                           timeout=timeout)]

    def cleanup(self, timeout=2000):
        """remove all FEMAG files in working directory
        (FEMAG 8.5 Rev 3282 or greater only)"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'cleanup'],
                                           timeout=timeout)]

    def release(self):
        """signal finish calculation task to load balancer to free resources
        (Docker Cloud environment only)
        """
        return [r.decode('latin1')
                for r in self.send_request(['close'], timeout=1000)]

    def info(self, timeout=2000):
        """get various resource information
        (FEMAG 8.5 Rev 3282 or greater only)"""
        response = [r.decode('latin1')
                    for r in self.send_request(['CONTROL', 'info'],
                                               timeout=timeout)]
        status = json.loads(response[0])['status']
        if status == 'ok':
            try:
                self.femaghost = json.loads(response[1])['addr']
                logger.info("Set femaghost %s", self.femaghost)
            except KeyError:
                # ignore femaghost setting if no addr is returned (Windows only)
                pass
        return response

    def publishLevel(self, level):
        """set publish level"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'publish = {}'
                                            .format(level)], timeout=10000)]

    def getfile(self, filename=''):
        """get file (FEMAG 8.5 Rev 3282 or greater only)"""
        response = self.send_request(
            ['CONTROL', f'getfile = {filename}'], timeout=1000)
        return [response[0].decode('latin1'),
                response[1] if len(response) else b'']

    def exportsvg(self, fslcmds, timeout=10000):
        """get svg format from fsl commands (if any graphic created)
        (since FEMAG 8.5 Rev 3343) """
        response = self.send_request(['SVG', fslcmds], timeout=timeout)
        rc = json.loads(response[0].decode('latin1'))
        if rc['status'] == 'ok':
            return self.getfile(rc['result_file'][0])
        return [s.decode('latin1') for s in response]

    def interrupt(self):
        """send push message to control port to stop current calculation"""
        context = zmq.Context.instance()
        if not self.femaghost:
            self.femaghost = '127.0.0.1'
        ctrl = context.socket(zmq.PUSH)
        ctrl.connect('tcp://{0}:{1}'.format(
            self.femaghost, self.port+2))
        logger.info("Interrupt %s", self.femaghost)
        ctrl.send_string('interrupt')
        ctrl.close()

    def copy_winding_file(self, name, wdg):
        wdg.write(name, self.workdir)
        self.upload(os.path.join(self.workdir, name+'.WID'))

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

    def __call__(self, machine, simulation):
        """setup fsl file, run calculation and return BCH results
        Args:
          machine: dict with machine parameters or name of model
          simulation; dict with simulation parameters

        Raises:
           FemagError
        """
        if isinstance(machine, str):
            modelpars = dict(name=machine)
        else:
            modelpars = machine
        if 'exit_on_end' not in modelpars:
            modelpars['exit_on_end'] = 'false'
        if 'exit_on_error' not in modelpars:
            modelpars['exit_on_error'] = 'false'
        response = self.send_fsl(['save_model("close")'] +
                                 self.create_fsl(modelpars,
                                                 simulation))
        r = json.loads(response[0])
        if r['status'] != 'ok':
            raise FemagError(r['message'])

        if not simulation:
            return [r, dict(name=machine['name'])]

        result_file = r['result_file'][0]
        if simulation['calculationMode'] == "pm_sym_loss":
            return self.read_los(self.modelname)

        status, content = self.getfile(result_file)
        r = json.loads(status)
        if r['status'] == 'ok':
            bch = femagtools.bch.Reader()
            bch.read(content.decode('latin1'))
            if simulation['calculationMode'] == 'pm_sym_fast':
                if simulation.get('shortCircuit', False):
                    logger.info("Short Circuit")
                    simulation.update(
                        get_shortCircuit_parameters(bch,
                                                    simulation.get('initial', 2)))
                    builder = femagtools.fsl.Builder()
                    response = self.send_fsl(
                        builder.create_shortcircuit(simulation))
                    r = json.loads(response[0])
                    if r['status'] != 'ok':
                        raise FemagError(r['message'])
                    result_file = r['result_file'][0]
                    status, content = self.getfile(result_file)
                    r = json.loads(status)
                    if r['status'] == 'ok':
                        bch.read(content.decode('latin1'))

            return bch
        raise FemagError(r['message'])
