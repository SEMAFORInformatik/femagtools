"""
    femagtools.docker
    ~~~~~~~~~~~~~~~~~

    Running FEMAG on Docker Swarm


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


def get_port_binding():
    """returns list with dict(HostIp, HostPort) of all
    running femag containers:
      [{'HostPort': '25555', 'HostIp': '0.0.0.0'}, 
       {'HostPort': '15555', 'HostIp': '0.0.0.0'}, 
       {'HostPort': '5555', 'HostIp': '0.0.0.0'}]
    """
    import docker
    client = docker.from_env()
    return [c.attrs['NetworkSettings']['Ports']['5555/tcp'][0]
            for c in client.containers.list(
                    filters={'label': 'org.label-schema.name=profemag/femag'})]


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
    def __init__(self, queue, workdir, port, host):
        threading.Thread.__init__(self)
        self.queue = queue
        self.container = femagtools.femag.ZmqFemag(
            workdir,
            port, host)
        
    def run(self):
        while True:
            task = self.queue.get()
            if task is None:
                break
            fslfile = os.path.join(task.directory, task.fsl_file)
            logger.info('Docker task %s %s', task.id, task.fsl_file)
            fslcmds = ['save_model(close)',
                       "chdir('{}')".format(
                           os.path.split(task.directory)[-1])]
            with open(fslfile) as f:
                fslcmds += f.readlines()
            r = [json.loads(s)
                 for s in self.container.send_fsl(
                         '\n'.join(fslcmds), publish_receive)[:2]]
            try:
                if r[0]['status'] == 'ok':
                    task.status = 'C'
                else:
                    task.status = 'X'
            except:
                task.status = 'X'
            logger.info("Finished %s", r)
            self.queue.task_done()
        
    
class Engine(object):

    """The Docker Engine

       execute Femag-Simulations with docker
    """
    def __init__(self, hosts=[]):
        self.hosts = hosts
        self.femag_port = int(os.environ.get('FEMAG_PORT', 5555))

    def create_job(self, workdir):
        """Create a FEMAG :py:class:`CloudJob`

        Args:
            workdir (str): The workdir where the calculation files are stored

        Return:
            job (:class:`Job`)
        """
        if not self.hosts:
            raise ValueError("empty host list")
        self.queue = Queue()
        self.async_femags = [AsyncFemag(self.queue,
                                        workdir,
                                        self.femag_port,
                                        h)
                             for h in self.hosts]
        
        self.job = femagtools.job.Job(workdir)
        return self.job

    def submit(self):
        """Starts the FEMAG calculation(s) as Docker containers

        Return:
            number of started tasks (int)
        """
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
        for async_femag in self.async_femags:
            async_femag.join()
            
        return [t.status for t in self.job.tasks]


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    engine = Engine()
    job = engine.create_job('/tmp/tar/docker-femag')
    for _ in range(3):
        t = job.add_task()
        t.add_file('femag.fsl',
                   content=[''])
    engine.submit()
    status = engine.join()
    print("Status {}".format(status))
    for t in engine.job.tasks:
        print(t.directory)
    
