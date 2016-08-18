"""
    femagtools.engin.multiproc
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    Creating and managing multicore/multiprocessing jobs

    :copyright: (c) (c) 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import multiprocessing
import subprocess
import os
import logging
from .job import Job

logger = logging.getLogger(__name__)


def run_femag(workdir, fslfile):
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


class MultiProc:
    def __init__(self, process_count=multiprocessing.cpu_count()):
        self.process_count = process_count
        
    def createJob(self, workdir):
        self.job = Job(workdir)
        return self.job
    
    def submit(self):
        pool = multiprocessing.Pool(self.process_count)
        self.tasks = [pool.apply_async(run_femag,
                                       args=(t.directory, t.fsl_file))
                      for t in self.job.tasks]
        return len(self.tasks)

    def join(self):
        results = [task.get() for task in self.tasks]
        status = []
        for t, r in zip(self.job.tasks, results):
            t.status = 'C' if r == 0 else 'X'
            status.append(t.status)
            
        return status

    def dumpData(self, workdir, parameters):
        '''dump data for later analysis'''
        pass
