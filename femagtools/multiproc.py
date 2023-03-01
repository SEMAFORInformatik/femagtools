"""
    femagtools.engine.multiproc
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Creating and managing multicore/multiprocessing jobs



    :authors: R. Tanner, N. Mauchle
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
try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

logger = logging.getLogger(__name__)

numpat = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)\s*')


class ProtFile:
    def __init__(self, dirname, num_cur_steps):
        self.size = 0
        self.looplen = 0
        self.cur_steps = [1, num_cur_steps]
        self.n = 0
        self.num_loops = 0
        self.dirname = dirname

    def percent(self):
        if self.looplen > 0:
            return 100 * self.n / self.looplen
        return 0

    def update(self):
        p = list(pathlib.Path(self.dirname).glob('*.PROT'))
        if p:
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
                return ''

        return f'{self.percent():3.1f}%'  # {self.n}/{self.looplen}'


class ProgressLogger(threading.Thread):
    def __init__(self, dirs, num_cur_steps, timestep):
        threading.Thread.__init__(self)
        self.dirs = dirs
        self.num_cur_steps = num_cur_steps
        self.running = False
        self.timestep = timestep

    def run(self):
        self.running = True
        protfiles = [ProtFile(d, self.num_cur_steps)
                     for d in self.dirs]
        while self.running:
            time.sleep(self.timestep)
            logmsg = [p.update() for p in protfiles]
            summary = [l  # f'<{i}> {l}'
                       for i, l in enumerate(logmsg)
                       if l]
            if summary:
                logger.info('Samples %s', ', '.join(summary))
            else:
                logger.info('collecting FE losses ...')
                return

    def stop(self):
        self.running = False


def run_femag(cmd, workdir, fslfile):
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
            proc = subprocess.Popen(cmd + ['-b', fslfile],
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

            logger.info("Finished pid: %d return %d",
                        proc.pid, proc.returncode)
            return proc.returncode
        except OSError as e:
            logger.error("Starting process failed: %s, Command: %s", e, cmd)
            raise

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
        progress_timestep: time step in seconds for progress log messages if > 0)
    """

    def __init__(self, **kwargs):
        self.process_count = kwargs.get('process_count', None)
        cmd = kwargs.get('cmd', '')
        if cmd:
            self.cmd = [cmd]
            if platform.system() == 'Windows':
                self.cmd.append('-m')

        self.progressLogger = 0
        self.progress_timestep = kwargs.get('timestep', 3)

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

        self.pool = multiprocessing.Pool(self.process_count)
        self.tasks = [self.pool.apply_async(run_femag,
                                            args=(t.cmd,
                                                  t.directory,
                                                  t.fsl_file))
                      for t in self.job.tasks]
        self.pool.close()

        if (self.progress_timestep and
                self.job.num_cur_steps):
            self.progressLogger = ProgressLogger(
                [t.directory for t in self.job.tasks],
                num_cur_steps=self.job.num_cur_steps,
                timestep=self.progress_timestep)
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
                    logger.error(t.errmsg)
            status.append(t.status)
        if self.progressLogger:
            self.progressLogger.stop()
        return status

    def terminate(self):
        logger.info("terminate Engine")
        if self.progressLogger:
            self.progressLogger.stop()
        # terminate pool
        try:
            self.pool.terminate()
            self.pool.close()
        except AttributeError:
            logger.warn("%s", e)

