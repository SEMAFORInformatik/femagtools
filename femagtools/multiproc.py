"""
    femagtools.engin.multiproc
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    Creating and managing multicore/multiprocessing jobs

    :copyright: (c) (c) 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
    :authors: R. Tanner, N. Mauchle
"""
import multiprocessing
import subprocess
import os
import logging
from .job import Job

logger = logging.getLogger(__name__)


def run_femag(workdir, fslfile):
    """Start the femag command as subprocess.

    :internal:

    :param workdir: The workdir where the calculation files are stored
    :param fslfile: The name of the start file (usually femag.fsl)
    """
    logger.info('FEMAG %s: %s', workdir, fslfile)
    with open(os.path.join(workdir, "femag.out"), "wb") as out, \
         open(os.path.join(workdir, "femag.err"), "wb") as err:
        proc = subprocess.Popen(['xfemag', '-b', fslfile],
                                shell=False,
                                stdout=out,
                                stderr=err,
                                cwd=workdir)
        # wait
        proc.wait()

    logger.info("Finished pid: %d return %d", proc.pid, proc.returncode)
    return proc.returncode


class Engine:
    """The MultiProc engine uses a pool of local calculation processes.

    This is more or less a decorator for the `Python multiprocessing Module
    <https://docs.python.org/3.6/library/multiprocessing.html>`_

    :param process_count: number of processes (cpu_count() if None)
    """
    def __init__(self, process_count=None):
        self.process_count = process_count

    def create_job(self, workdir):
        """Create a FEMAG :py:class:`Job`

        :param workdir: The workdir where the calculation files are stored
        :return: FEMAG :py:class:`Job`
        """
        self.job = Job(workdir)
        return self.job

    def submit(self):
        """Starts the FEMAG calculation(s) with the internal
        :py:meth:`multiproc.run_femag` function

        :return: length of started tasks
        """
        pool = multiprocessing.Pool(self.process_count)
        self.tasks = [pool.apply_async(run_femag,
                                       args=(t.directory, t.fsl_file))
                      for t in self.job.tasks]
        return len(self.tasks)

    def join(self):
        """Wait until all calculations are finished

        :return: list of all calculations status (C = Ok, X = error)
        """
        results = [task.get() for task in self.tasks]
        status = []
        for t, r in zip(self.job.tasks, results):
            t.status = 'C' if r == 0 else 'X'
            status.append(t.status)

        return status
