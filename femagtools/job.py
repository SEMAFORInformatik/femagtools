# -*- coding: utf-8 -*-
"""
    femagtools.job
    ~~~~~~~~~~~~~~

    Creating and managing calculation jobs.



"""
import os
import platform
import shutil
import glob
import femagtools.bch
import femagtools.config as cfg
import logging
import uuid

logger = logging.getLogger(__name__)


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
        """adds a file required by this task

        Args:
            fname: file name 
            content: list of str written to file if not None"""
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
        result = femagtools.bch.Reader()
        # read latest bch file if any
        bchfile_list = sorted(glob.glob(os.path.join(
            self.directory, '*_[0-9][0-9][0-9].B*CH')))
        if bchfile_list:
            with open(bchfile_list[-1]) as f:
                logger.info("Reading %s",
                            bchfile_list[-1])
                result.read(f)
            #logger.info("%s %s", result.version, result.type),
        else:
            msg = 'no BCH files in {}'.format(self.directory)
            logger.error(msg)
            result = dict(error=msg)
        return result

    def __eq__(self, other):
        return isinstance(other, type(self)) and \
            self.directory == other.directory


class CloudTask(Task):
    def __init__(self, id, directory):
        super(self.__class__, self).__init__(id, directory)
        import tarfile
        self.file = "{}.tar.gz".format(self.directory)
        self.tar_file = tarfile.open(self.file, "w:gz")
        # Used for amazon
        self.ec2_instance = None

    def add_file(self, fname, content=None):
        base = os.path.basename(fname)
        self.transfer_files.append(fname)
        if os.path.splitext(base)[-1] == '.fsl':
            self.fsl_file = base

        info = self.tar_file.tarinfo()
        info.name = base
        if content is None:
            info.size = os.path.getsize(fname)
            self.tar_file.addfile(info, open(fname, 'rb'))
            return

        import io
        if type(content) is list:
            content = '\n'.join(content)
        info.size = len(content)
        data = io.BytesIO(str.encode(content))
        self.tar_file.addfile(info, data)


class Job(object):
    """represents a FEMAG job consisting of one or more tasks
    each to be executed by a dedicated process
    """
    def __init__(self, basedir):
        self.runDirPrefix = ''
        self.basedir = basedir
        self.tasks = []

    def cleanup(self):
        """removes all task directories of previous run"""
        for task in self.tasks:
            logger.info("rm %s", task.directory)
            shutil.rmtree(task.directory, ignore_errors=True)
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

        OpSys = "WINDOWS" if platform.system() == "Windows" else "LINUX"

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
            'Executable   = {}'.format(cfg.get_executable()),
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


class CloudJob(Job):
    """Inheritance of :py:class:`Job`

    Represents a femag amazon job"""
    def __init__(self, basedir):
        super(self.__class__, self).__init__(basedir)

    def add_task(self):
        "adds a new :py:class:`AmazonTask` to this job"
        taskid = "{}-{}".format(str(uuid.uuid4()), len(self.tasks))
        dir = os.path.join(self.basedir,
                           '{}{:d}'.format(self.runDirPrefix,
                                           len(self.tasks)))
        self.tasks.append(CloudTask(taskid, dir))
        return self.tasks[-1]
