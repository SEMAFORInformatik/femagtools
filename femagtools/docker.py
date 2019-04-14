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
    def __init__(self, queue, port, host):
        threading.Thread.__init__(self)
        self.queue = queue
        self.container = femagtools.femag.ZmqFemag(
            port, host)
        
    def run(self):
        """execute femag fsl task in task directory"""
        while True:
            task = self.queue.get()
            if task is None:
                break

            r = self.container.cleanup()
            for f in task.transfer_files:
                if f != task.fsl_file:
                    r = self.container.upload(
                        os.path.join(task.directory, f))
            fslfile = os.path.join(task.directory, task.fsl_file)
            logger.info('Docker task %s %s',
                        task.id, task.fsl_file)
            fslcmds = []
            with open(fslfile) as f:
                fslcmds = f.readlines()
            ret = self.container.send_fsl(fslcmds +
                                          ['save_model(close)'])
            # TODO: add publish_receive
            logger.debug(ret)
            r = [json.loads(s) for s in ret]
            logger.debug("Finished %s", r)
            try:
                if r[0]['status'] == 'ok':
                    task.status = 'C'
                    bchfile = r[0]['result_file'][0]
                    status, content = self.container.getfile(bchfile)
                    logging.info("get results %s: status %s len %d",
                                 task.id, status, len(content))
                    with open(os.path.join(task.directory,
                                           bchfile), 'wb') as f:
                        f.write(content)
                else:
                    task.status = 'X'
                    logger.warn("%s: %s", task.id, r['message'])
            except (KeyError, IndexError):
                task.status = 'X'

            logger.info("Task %s end status %s",
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

    def submit(self):
        """Starts the FEMAG calculation(s) as Docker containers

        Return:
            number of started tasks (int)
        """
        self.queue = Queue()
        self.async_femags = [AsyncFemag(self.queue,
                                        self.port, self.dispatcher)
                             for i in range(self.num_threads)]

        for async_femag in self.async_femags:
            async_femag.start()

        for task in self.job.tasks:
            self.queue.put(task)
            
        return len(self.job.tasks)

    def join(self):
        """Wait until all calculations are finished

        Return:
            list of all calculations status (C = Ok, X = error) (:obj:`list`)
        """
        # block until all tasks are done
        self.queue.join()

        # stop workers
        for _ in self.async_femags:
            self.queue.put(None)

        # wait for all workers
        for async_femag in self.async_femags:
            async_femag.join()

        return [t.status for t in self.job.tasks]
