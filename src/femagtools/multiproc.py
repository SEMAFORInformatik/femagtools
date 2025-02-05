"""manage multicore/multiprocessing jobs

"""
import platform
import multiprocessing
import subprocess
import time
import os
import re
import threading
import pathlib
import logging
from .job import Job
import femagtools.config as cfg
import femagtools.zmq
from femagtools.zmq import SubscriberTask
try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

logger = logging.getLogger(__name__)

class LicenseError(Exception):
    pass

class ProgressLogger(threading.Thread):
    def __init__(self, dirs, num_cur_steps, timestep, notify):
        threading.Thread.__init__(self)
        self.protfiles = [femagtools.zmq.ProtFile(d, num_cur_steps)
                          for d in dirs]
        self.numTot = len(dirs)
        self.running = False
        self.timestep = timestep
        self.notify = notify

    def run(self):
        self.running = True
        while self.running:
            if self.timestep > 0:
                time.sleep(self.timestep)
            logmsg = [p.update() for p in self.protfiles]
            summary = [l  # f'<{i}> {l}'
                       for i, l in enumerate(logmsg)
                       if l]
            if summary:
                labels = set([p.name for p in self.protfiles])
                logger.info('%s: %s',
                            ', '.join(labels),
                            ', '.join(summary))
                if self.notify:
                    numOf = f"{summary.count('100.0%')} of {self.numTot}"
                    percent = sum([float(i[:-1])
                                   for i in summary]) / self.numTot
                    self.notify(
                        ["progress_logger",
                         f"{self.numTot}:{numOf}:{percent}:{' '.join(summary)}"])
            else:
                # TODO: log message might be misleading
                logger.debug('collecting FE losses ...')
                return

    def stop(self):
        self.running = False


def run_femag(cmd, workdir, fslfile, port):
    """Start the femag command as subprocess.

    :internal:

    Args:
        cmd: (list) The program (executable image) to be run
        workdir: The workdir where the calculation files are stored
        fslfile: The name of the start file (usually femag.fsl)
    """
    with open(os.path.join(workdir, "femag.out"), "wb") as out, \
            open(os.path.join(workdir, "femag.err"), "wb") as err:
        try:
            args = ['-b', str(port), fslfile] if port else ['-b', fslfile]
            proc = subprocess.Popen(cmd + args,
                                    shell=False,
                                    stdin=DEVNULL,
                                    stdout=out,
                                    stderr=err,
                                    cwd=workdir)
            logger.info('%s (pid %d, workdir %s)', cmd, proc.pid, workdir)
            # write pid file
            with open(os.path.join(workdir, 'femag.pid'), 'w') as pidfile:
                pidfile.write("{}\n".format(proc.pid))

            # wait
            proc.wait()
            os.remove(os.path.join(workdir, 'femag.pid'))

            #logger.info("Finished pid: %d return %d",
            #            proc.pid, proc.returncode)
            #return proc.returncode
        except OSError as e:
            logger.error("Starting process failed: %s, Command: %s", e, cmd)
            raise

    # raise License Error
    if proc.returncode != 0:
        with open(os.path.join(workdir, "femag.err"), "r") as err:
            for line in err:
                if 'license' in line:
                    raise LicenseError(line)

    logger.info("Finished pid: %d return %d", proc.pid, proc.returncode)
    return proc.returncode


class Engine:
    """The MultiProc engine uses a pool of local calculation processes.

    This is more or less a decorator for the `Python multiprocessing Module
    <https://docs.python.org/3.6/library/multiprocessing.html>`_

    Args:
        cmd: the program (executable image) to be run
            (femag dc is used if None)
        process_count: number of processes (cpu_count() if None)
        timestep: time step in seconds for progress log messages if > 0)
    """

    def __init__(self, **kwargs):
        self.process_count = kwargs.get('process_count', None)
        self.notify = kwargs.get('notify', None)
        # cogg_calc mode, subscribe xyplot
        self.calc_mode = kwargs.get('calc_mode')
        self.port = kwargs.get('port', 0)
        self.curve_label = kwargs.get('curve_label')
        cmd = kwargs.get('cmd', '')
        if cmd:
            self.cmd = [cmd]
            if platform.system() == 'Windows':
                self.cmd.append('-m')

        self.progressLogger = 0
        self.progress_timestep = kwargs.get('timestep', -1)
        self.subscriber = None
        self.job = None
        self.tasks = []

    def create_job(self, workdir):
        """Create a FEMAG :py:class:`Job`

        Args:
            workdir: The workdir where the calculation files are stored

        Return:
            FEMAG :py:class:`Job`
        """
        self.job = Job(workdir)
        return self.job

    def submit(self, extra_result_files=[]):
        """Starts the FEMAG calculation(s) with the internal
        :py:meth:`multiproc.run_femag` function

        Return:
            length of started tasks
        """
        # must check if cmd is set:
        args = []
        if platform.system() == 'Windows':
            args.append('-m')

        try:
            for t in self.job.tasks:
                t.cmd = self.cmd
        except AttributeError:
            for t in self.job.tasks:
                t.cmd = [cfg.get_executable(
                    t.stateofproblem)] + args

        num_proc = self.process_count
        if not num_proc and multiprocessing.cpu_count() > 1:
            num_proc = min(multiprocessing.cpu_count()-1, len(self.job.tasks))
        self.pool = multiprocessing.Pool(num_proc)
        if self.port:
            header = [b'progress']
            if self.calc_mode == 'cogg_calc':
                header +=[b'xyplot']
            self.subscriber = [SubscriberTask(port=self.port + i * 5,
                                              host='127.0.0.1',
                                              notify=self.notify,
                                              header=header,
                                              curve_label=self.curve_label,
                                              num_cur_steps=self.job.num_cur_steps,
                                              timestep=self.progress_timestep
                                              )
                               for i, t in enumerate(self.job.tasks)]
            [s.start() for s in self.subscriber]
            self.tasks = [self.pool.apply_async(
                run_femag, args=(t.cmd, t.directory, t.fsl_file, self.port + i * 5))
                          for i, t in enumerate(self.job.tasks)]
        else:
            self.tasks = [self.pool.apply_async(
                run_femag, args=(t.cmd, t.directory, t.fsl_file, 0))
                          for t in self.job.tasks]
        self.pool.close()

        # only works on linux
        if (self.progress_timestep and not self.port and
            self.job.num_cur_steps):
            self.progressLogger = ProgressLogger(
                [t.directory for t in self.job.tasks],
                num_cur_steps=self.job.num_cur_steps,
                timestep=self.progress_timestep,
                notify=self.notify)
            self.progressLogger.start()
        return len(self.tasks)

    def join(self):
        """Wait until all calculations are finished

        Return:
            list of all calculations status (C = Ok, X = error)
        """
        exitcodes = [task.get() for task in self.tasks]
        status = []
        for t, ec in zip(self.job.tasks, exitcodes):
            t.status = 'C'
            if ec != 0:
                t.status = 'X'
                errmsg = pathlib.Path(t.directory) / 'femag.err'
                if errmsg.exists():
                    t.errmsg = errmsg.read_text()
                    if t.errmsg:
                        logger.error(t.errmsg)
            status.append(t.status)
        self.stopThreads()
        self.pool = None # garbage collector deletes threads
        return status

    def stopThreads(self):
        """ stop all running treads
        """
        if self.progressLogger:
            self.progressLogger.stop()
        if self.port and self.subscriber:
            [s.stop() for s in self.subscriber]
            SubscriberTask.clear()
            self.subscriber = None

    def terminate(self):
        """ terminate all
        """
        logger.info("terminate Engine")
        self.stopThreads()

        # terminate pool
        try:
            if self.pool:
                self.pool.terminate()
                self.pool = None # garbage collector deletes threads
        except AttributeError as e:
            logger.warn("%s", e)
