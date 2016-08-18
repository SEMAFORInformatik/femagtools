# -*- coding: utf-8 -*-
"""
    femagtools.job
    ~~~~~~~~~~~~~~

    Creating and managing calculation jobs.

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import os
import platform
import shutil
import glob
import femagtools.bch
import logging
import uuid
try:
    import boto3
except:
    pass

logger = logging.getLogger(__name__)


class JobError(Exception):
    """raised when an error occured"""


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


class Task(object):
    """represents a single execution unit that may include data files"""
    def __init__(self, id, directory):
        self.directory = directory
        try:
            os.makedirs(self.directory)
        except OSError as e:
            if e.errno != 17:
                e.args += self.directory
                raise
            pass
        self.transfer_files = []
        self.status = None
        self.fsl_file = None
        self.id = id
        
    def add_file(self, fname, content=None):
        """adds a file required by this task"""
        base = os.path.basename(fname)
        self.transfer_files.append(fname)
        if os.path.splitext(base)[-1] == '.fsl':
            self.fsl_file = base

        if content is None:
            dest = os.path.join(self.directory, base)
            if not os.access(dest, os.R_OK):
                shutil.copy(fname, dest)
            return

        # this file has to be created
        with open(os.path.join(self.directory, base), 'w') as f:
            f.writelines('\n'.join(content))

    def get_results(self):
        """returns result of most recent BCH file"""
        result = femagtools.bch.BchReader()
        # read latest bch file if any
        bchfile_list = sorted(glob.glob(os.path.join(
            self.directory, '*_[0-9][0-9][0-9].B*CH')))
        if bchfile_list:
            with open(bchfile_list[-1]) as f:
                result.read(f)
            logger.info("%s %s", result.version, result.type),
        else:
            msg = 'no BCH files in {}'.format(self.directory)
            logger.error(msg)
            result = dict(error=msg)
        return result

    def __eq__(self, other):
        return isinstance(other, type(self)) and \
            self.directory == other.directory


class AmazonTask(Task):
    def __init__(self, id, directory):
        super(self.__class__, self).__init__(id, directory)
        # These attributes are set during runtime
        self.exit_code = None
        self.ec2_instance = None
        self.container_id = None
        self.amazon_task_id = None
        self.sync_command_id = None
        self.finished = False

        def __str__(self):
            return ", ".join(('id: {}',
                              'directory: {}',
                              'ec2_instance: {}',
                              'exit_code: {}',
                              'container_id {}',
                              'task_id: {}',
                              'sync_command_id: {}',
                              'finished: {}')).format(
                self.id,
                self.directory,
                self.ec2_instance,
                self.exit_code,
                self.container_id,
                self.task_id,
                self.sync_command_id,
                self.finished)

        def __repr__(self):
            return self.__str__()


class Job(object):
    """represents a FEMAG job consisting of one or more tasks
    each to be executed by a dedicated process
    """
    def __init__(self, basedir):
        self.runDirPrefix = ''
        self.basedir = basedir
        self.tasks = []
        
    def cleanup(self):
        """removes all files and directories of previous run"""
        try:
            os.makedirs(self.basedir)
        except:
            if not os.path.isdir(self.basedir):
                e.args += self.basedir
                raise
        for d in os.listdir(self.basedir):
            d = os.path.join(self.basedir, d)
            if os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)
        self.tasks = []
        
    def add_task(self):
        "adds a new task to this job"
        taskid = "{}-{}".format(str(uuid.uuid4()), len(self.tasks))
        dir = os.path.join(self.basedir,
                           '{}{:d}'.format(self.runDirPrefix,
                                           len(self.tasks)))
        self.tasks.append(Task(taskid, dir))
        return self.tasks[-1]
    
    def setExitStatus(self, taskid, status):
        "set exit status of task"
        self.tasks[taskid].status = status

    def get_results(self):
        for t in self.tasks:
            yield t.get_results()


class CondorJob(Job):
    """represents a femag job that is to be run in HT Condor"""
    def __init__(self, basedir):
        super(self.__class__, self).__init__(basedir)
    
    def prepareDescription(self):
        # create a flatten list of all files to be transferred
        transfer_files = [item for sublist in
                          [t.transfer_files for t in self.tasks]
                          for item in sublist]
        if platform.system() == "Windows":
            wfemag = which("wfemagw64.exe")
            if not wfemag:
                raise JobError("wfemag not found in {}".format(
                    os.environ["PATH"]))
            for l in ['libstdc++-6.dll', 'libgcc_s_dw2-1.dll']:
                lib = which(l)
                if not lib:
                    raise JobError("{} not found in {}".format(
                        l, os.environ["PATH"]))
                
                transfer_files.append(lib)
            Executable = wfemag
            OpSys = "WINDOWS"

        else:
            xfemag = which("xfemag")
            if not xfemag:
                raise JobError("xfemag not found in {}".format(
                    os.environ["PATH"]))
            Executable = xfemag
            OpSys = "LINUX"

        fslFilename = ''
        for f in transfer_files:
            b, ext = os.path.splitext(f)
            if ext == '.fsl':
                fslFilename = f
                break
        numRuns = len(self.tasks)
# $(Process)                os.path.join(self.basedir, self.runDirPrefix)),
        submit = [
            'InitialDir   = {}/{}$(Process)'.format(
                self.basedir, self.runDirPrefix),
            'Universe     = vanilla',
            'Executable   = {}'.format(Executable),
            'Arguments    = -b {}'.format(fslFilename),
            'Output       = femag.out',
            'Log          = femag.log',
            'Error        = femag.err',
            'Notification = never',
            'input        = /dev/null',
            'transfer_input_files = {}'.format(','.join(set(transfer_files))),
            'requirements = (OpSys=="{}") && Arch=="x86_64"'.format(OpSys),
            'Queue {}'.format(numRuns), '']
        filename = os.path.join(self.basedir, "femag.submit")
        with open(os.path.join(self.basedir,
                               "femag.submit"), 'w') as submitFile:
            submitFile.writelines('\n'.join(submit))
        return filename


class AmazonJob(Job):
    """represents a femag job that is to be run in HT Condor"""
    def __init__(self, basedir):
        super(self.__class__, self).__init__(basedir)

    def cleanup(self):
        super(self.__class__, self).cleanup()
        # Amazon cleanup
        s3 = boto3.resource('s3')
        for task in self.tasks:
            s3.buckets(task.id).objects.delete()
            s3.meta.client.delete_bucket(Bucket=task.id)

    def add_task(self):
        "adds a new task to this job"
        taskid = "{}-{}".format(str(uuid.uuid4()), len(self.tasks))
        dir = os.path.join(self.basedir,
                           '{}{:d}'.format(self.runDirPrefix,
                                           len(self.tasks)))
        self.tasks.append(AmazonTask(taskid, dir))
        return self.tasks[-1]


if __name__ == "__main__":
    import tempfile
    workdir = tempfile.mkdtemp()
    job = CondorJob(workdir)
    task = job.add_task()
    task.add_file('femag.fsl', ['exit_on_end=True'])
    job.prepareDescription()
    print("Done: {}".format(workdir))
