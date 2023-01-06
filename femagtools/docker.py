"""
    femagtools.docker
    ~~~~~~~~~~~~~~~~~

    Running FEMAG on Docker/Cloud


"""
import os
import json
import logging
import threading
import femagtools.femag
import femagtools.job
import time
try:
    from queue import Queue
except ImportError:
    from Queue import Queue  # python 2.7


logger = logging.getLogger(__name__)


def publish_receive(message):
    """handle messages from femag publisher"""
    topic, content = message  # "femag_log" + text
    # topics: femag_log, progress, file_modified,
    #   model_image, calc_image, field_image, babs_image, demag_image, color_scale
    if topic == 'femag_log' or topic == 'progress':
        logger.info("%s: %s", topic, content.strip())
    else:
        logger.info('%s: len %d', topic, len(content.strip()))


class AsyncFemag(threading.Thread):
    def __init__(self, queue, port, host,
                 extra_result_files):
        threading.Thread.__init__(self)
        self.queue = queue
        self.container = femagtools.femag.ZmqFemag(
            port, host)
        self.extra_result_files = extra_result_files

    def _do_task(self, task):
        resend = True
        num_tries = 0
        while resend:
            r = self.container.cleanup(timeout=None)
            status = json.loads(r[0])
            resend = status['status'] == 'resend'
            if resend:
                time.sleep(1)
            num_tries += 1
            if num_tries > 4:
                logger.warning(
                    'Docker task %s cleanup status tries %d',
                    task.id, num_tries)
                return [dict(status='error', msg='too many resend')]

        if status['status'] != 'ok':
            logger.warning('Docker task %s cleanup status %s',
                           task.id, status)
            return [status]
        for f in task.transfer_files:
            if f != task.fsl_file:
                r = self.container.upload(
                    os.path.join(task.directory,
                                 os.path.basename(f)))
                status = json.loads(r[0])
                if status['status'] != 'ok':
                    logger.warning('Docker task %s upload status %s',
                                   task.id, status)
                    return [status]
            fslfile = os.path.join(task.directory, task.fsl_file)
        fslcmds = []
        with open(fslfile) as f:
            fslcmds = f.readlines()
        ret = self.container.send_fsl(fslcmds +
                                      ['save_model(close)'])
        # TODO: add publish_receive
        try:
            return [json.loads(s) for s in ret]
        except Exception as e:
            logger.error("%s: %s", e, ret, exc_current=True)
            return [{'status': 'error', 'message': str(ret)}]

    def run(self):
        """execute femag fsl task in task directory"""
        while True:
            task = self.queue.get()
            if task is None:
                break

            try:
                r = self._do_task(task)
                if r[0]['status'] == 'ok':
                    task.status = 'C'
                    result_files = []
                    if 'result_file' in r[0]:
                        result_files = [r[0]['result_file'][0]]
                    result_files += self.extra_result_files + task.extra_result_files
                    for fname in result_files:
                        status, content = self.container.getfile(fname)
                        logging.debug("get results %s: status %s len %d",
                                      task.id, status, len(content))
                        with open(os.path.join(task.directory,
                                               fname), 'wb') as f:
                            f.write(content)
                else:
                    task.status = 'X'
                    logger.warning("%s: %s", task.id, r[0])
            except (KeyError, IndexError):
                task.status = 'X'
                logger.error("AsyncFemag", exc_info=True)

            logger.debug("Task %s end status %s",
                         task.id, task.status)
            ret = self.container.release()
            self.queue.task_done()
        self.container.close()


class Engine(object):

    """The Docker Engine

       execute Femag-Simulations with docker

       Args:
         dispatcher (str): hostname of dispatcher
         port (int): port number of dispatcher
         num_threads: number of threads to send requests
    """

    def __init__(self, dispatcher='127.0.0.1', port=5000,
                 num_threads=5):
        self.port = port
        self.dispatcher = dispatcher
        self.num_threads = num_threads
        self.async_femags = None

    def create_job(self, workdir):
        """Create a FEMAG :py:class:`CloudJob`

        Args:
            workdir (str): The workdir where the calculation files are stored

        Return:
            job (:class:`Job`)
        """
        self.job = femagtools.job.Job(workdir)
        return self.job

    def submit(self, extra_result_files=[]):
        """Starts the FEMAG calculation(s) as Docker containers

        Return:
            number of started tasks (int)
        """
        self.queue = Queue()
        for task in self.job.tasks:
            self.queue.put(task)

        logger.info("Request %d workers on %s:%d (num tasks %d)",
                    self.num_threads, self.dispatcher, self.port,
                    len(self.job.tasks))
        self.async_femags = [
            AsyncFemag(
                self.queue,
                self.port, self.dispatcher,
                extra_result_files)
            for i in range(self.num_threads)]

        for async_femag in self.async_femags:
            async_femag.start()

        return len(self.job.tasks)

    def join(self):
        """Wait until all calculations are finished

        Return:
            list of all calculations status (C = Ok, X = error) (:obj:`list`)
        """
        # block until all tasks are done
        logger.debug("join: block until all tasks are done")
        self.queue.join()

        # stop workers
        logger.debug("join: stop workers")
        for _ in self.async_femags:
            self.queue.put(None)

        # wait for all workers
        logger.debug("join: join workers")
        for async_femag in self.async_femags:
            async_femag.join()

        logger.debug("join: done")
        return [t.status for t in self.job.tasks]
