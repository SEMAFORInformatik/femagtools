"""
Author: Nicolas Mauchle <mauchle@semafor.ch>
To use this class you have to install:

- gcloud
- googleapiclient

Also you have to configure your google account with
  gcloud init

The google cloud module for femag in short:
1. Create for every calculation a bucket
2. Upload the files for the calculation in the specific bucket
3. Start for every calculation one instance with a femag image
4. Wait until a exit_code file is stored in the bucket
5. Terminate instance
6. Sync back the results files

The connection between the buckets and the calculation are stored in the job dict

"""
import logging
import os
import threading      # Not supported from httplib
import time
import femagtools
import random         # later used in create project
import pdb
from .config import Config

logger = logging.getLogger(__name__)

GOOGLE_MODULE_FOUND = True
try:
    from oauth2client.client import GoogleCredentials
    from googleapiclient import discovery
    from gcloud import storage
    from gcloud import resource_manager
except:
    GOOGLE_MODULE_FOUND = False


class IllegalNumberOfBuckets(Exception):
    def __init__(self):
        Exception.__init__(self)
        logger.error("Max 10 buckets are allowed")

class Engine():

    config_class = Config

    default_config = {
        'ENGINE': 'google',
        'SERVER_LOCATION': 'us-east1-b',  # 'f1-micro'
        'INSTANCE_TYPE': 'n1-standard-1',
        'IMAGE_ID': 'femag-v1',
        'FINISH_TASK_FILENAME': 'exit_code',
        'COMPANY_NAME': 'femag',
    }

    def __init__(self, buckets=None):
        """Initalize the Google cloud service

        Args:
            buckets:     to prevent uploads and use these buckets with data
        """
        if not GOOGLE_MODULE_FOUND:
            logger.error("Could not find one of [google-api-python-client, gcloud, googleapis-common-protos] modules")
            return

        if buckets and len(buckets) > 10:
            raise IllegalNumberOfBuckets()

        self.buckets = buckets
        self.job = None

        self.config = Config(self.default_config)

        # Create a new project
        self.project = self._create_new_project('ancient-ship-136307')

    def _create_data_buckets(self):
        """Every calculation has its own bucket with all data
        All buckets are created before any calculation.
        The name of a bucket is unique uuid number

        .. note::  We can only create 10 buckets in a short time period
        """

        # If we have already some buckets we have to assign
        # them to our jobs
        if self.buckets:
            logger.info("Buckets already exists")
            for idx, bucket in enumerate(self.buckets):
                self.job.tasks[idx].id = bucket['id']
                self.job.tasks[idx].directory = bucket['folder']
            return

        gcs = storage.Client(self.project.project_id)
        # Create new buckets for calculation
        logger.info("Create {} new buckets".format(len(self.job.tasks)))
        for t in self.job.tasks:
            gcs.create_bucket(t.id)

        logger.info("Buckets created")

    def _upload_files_to_buckets(self):
        """For every calculation there are some files to upload in the specific bucket
        To save time the upload is done in a seperate thread.

        .. note::

        Unfortunately the google storage model can not handle threads (SSL version error)
        To avoid that, we have to create a new storage client in every thread.
        """

        # Do not upload files if whe already have the buckets
        if self.buckets:
            logger.info("Files are already uploaded")
            return

        threads = []
        for t in self.job.tasks:
            logger.info("Start upload folder {} from task {}".format(t.file, t.id))
            thread = threading.Thread(target=self._upload, args=(t, ))
            thread.start()
            threads.append(thread)

        # Wait for all finished uploads
        self._wait_for_threads_finished(threads, "Uploading files")

    def _upload(self, task):
        """Upload the file to the google storage

        Args:
            task:  The task which define which folder (task.directory) should be uploaded
                      to which storage (task.id)
        """
        # We can not use the self.gcs cause we use it in a thread and if we use self.gcs an exception
        # is thrown -> Wrong SSL Version.
        gcs = storage.Client(self.project.project_id)

        # task.id is the uuid generated name for the bucket
        bucket = gcs.get_bucket(task.id)
        # Close the tar file
        task.tar_file.close()
        
        blob = storage.Blob(os.path.basename(task.file), bucket)
        with open(task.file, 'rb') as file:
                blob.upload_from_file(file)

        # Other possibility:
        # bucket = gcs.create_bucket('bucket-name')
        # blob = bucket.blob(filename)
        # blob.upload_from_string(uploaded_file.read(), content_type=uploaded_file.content_type )

    def _create_new_project(self, name):
        """At the moment it is not possible to create a new project overt the python api
        So we use always the default project

        Args:
            name (str):     The name of the project

        Return:
            Google Cloud project (:obj:`Google Cloud Project`)
        """
        gcrm = resource_manager.Client()

        # Create new project
        # project_id = "{}-{}".format(self.company, random.randint(1, 10000))
        # self.project = self.gcrm.new_project(project_id, name=name)
#        project.labels = {'tasks': ", ".join([t.id for t in self.job.tasks]), 'status': 'initalize'}
        # self.project.create()  # Create the project

        project = gcrm.fetch_project(name)
        logger.info("Created new project. ID: {}, Name: {}".format(project.project_id, project.name))
        return project

    def _start_instances(self):
        """ Prepare and start all instances for the calculation.

        """

        # Configuration for our machine. Here we put all our predefined configuration
        import json
        config = json.load(open(self.config['SERVER_CONFIG']))

        # If we want to use a debian image file from google:
        # image_response = self.compute.images().getFromFamily(project='debian-cloud',
        #                                                     family=self.image_id).execute()
        # source_disk_image = image_response['selfLink']

        # Use or own image disk with preinstalled femag and dependencies:
        for disk in config['disks']:
            if 'initializeParams' not in disk:
                disk['initializeParams'] = {'sourceImage': "projects/{}/global/images/{}".format(self.project.project_id, self.config['IMAGE_ID'])}

        if 'machineType' not in config:
            config['machineType'] = "zones/{}/machineTypes/{}".format(self.config['SERVER_LOCATION'],
                                                                      self.config['INSTANCE_TYPE'])

        config['metadata']['items'].append({'key': 'startup-script',
                                            'value': open(self.config['STARTUP_SCRIPT'], 'r').read()})
        # Startup instances and the calculation
        operations = []

        # You always have to be logged in. With GoogleCredentials.get_application_default()
        # it looks in the home directory for the credentials and use the default profile
        credentials = GoogleCredentials.get_application_default()

        # Get the compute object
        compute = discovery.build('compute', 'v1', credentials=credentials)
        for t in self.job.tasks:
            # Set same name for the instance as for the bucket.
            # Googlecloud does not allow instance name which starts with a number
            config['name'] = "{}-{}".format(self.config['COMPANY_NAME'].lower(), t.id)

            # Check if we defined the key 'bucket' in items
            # in loop 1 the bucket is not defined -> so we have to add it
            # after that we overwrite the bucket value
            # Every instance does download his bucket
            idx = [i for i, item in enumerate(config['metadata']['items']) if item['key'] == 'bucket']
            if idx:
                config['metadata']['items'][idx[0]] = {'key': 'bucket', 'value': t.id}
            else:
                config['metadata']['items'].append({'key': 'bucket', 'value': t.id})

            # Start the instance with our project_id zone and configuration
            operation = compute.instances().insert(
                project=self.project.project_id,
                zone=self.config['SERVER_LOCATION'],
                body=config).execute()

            operations.append(operation)

        # Wait until all instances are started up
        self._wait_for_operations_finished(operations)

    def _wait_for_threads_finished(self, threads, operation):
        """This generic methods waits until all threads are finished

        Args:
            list threads:    List of threads to check if they are finished
            str operation:   Name of the operation to write a meaningful log message

        """
        # Wait until all threads are not running
        while not all([t.isAlive() is False for t in threads]):
            time.sleep(5)

        # timer.cancel
        logger.info("{} are finished".format(operation))

    def _wait_for_operations_finished(self, operations):
        """This methods waits until all operations (started on google cloud) are finished

        Args:
            operations (:obj:`list`): A list of google cloud operations
        """
        threads = []
        # timeout = 300
        # timer = threading.Timer(timeout, thread.interrupt_main)
        for operation in [o['name'] for o in operations]:
            thread = threading.Thread(target=self._wait_for_operation, args=(operation, ))
            threads.append(thread)
            thread.start()

        # timer.start()
        self._wait_for_threads_finished(threads, "Operations")

    def _wait_for_operation(self, operation):
        """Wait until a operation is finished
        Is called from :py:meth:`_wait_for_operations_finished`

        Args:
            Operation One Google Operation Object
        """

        #  Do not use self.compute: SSL VERSION Exception when called in thread
        credentials = GoogleCredentials.get_application_default()
        compute = discovery.build('compute', 'v1', credentials=credentials)
        logger.info('Waiting for operation to finish...')
        while True:
            result = compute.zoneOperations().get(
                project=self.project.project_id,
                zone=self.config['SERVER_LOCATION'],
                operation=operation).execute()

            # Check if operation is DONE or has any error
            if result['status'] == 'DONE':
                if 'error' in result:
                    raise Exception(result['error'])

                return result

            time.sleep(10)

    def _delete_instances(self, instances):
        """Delete one ore more instances

        Args:
            instances (:obj:`list`):   A list of instances/task to be deleted

        """
        # You always have to be logged in. With GoogleCredentials.get_application_default()
        # it looks in the home directory for the credentials and use the default profile
        credentials = GoogleCredentials.get_application_default()

        # Get the compute object
        compute = discovery.build('compute', 'v1', credentials=credentials)
        operations = []
        for task in instances:
            logger.info("Delete instace: {}".format(task.id))
            # Be sure to make all strings lower cause of the google regex
            # (?:[a-z](?:[-a-z0-9]{0,61}[a-z0-9])?)
            instance_name = "{}-{}".format(self.config['COMPANY_NAME'].lower(), task.id)
            operations.append(compute.instances().delete(
                project=self.project.project_id,
                zone=self.config['SERVER_LOCATION'],
                instance=instance_name).execute())

        self._wait_for_operations_finished(operations)

    def _download_bucket(self, task):
        """Download all files from one bucket bucket name is task.id in a thread

        Args:
            task (:py:class:`CloudTask`):  The task which is assigned to this bucket to get the bucket name
        """
        # !Threads!
        gcs = storage.Client(self.project.project_id)
        bucket = gcs.get_bucket(task.id)

        if not os.path.exists(task.directory):
            os.makedirs(task.directory)

        # Get all files which are stored in our bucket
        files = bucket.list_blobs()
        for file in files:
            dest = "{}/{}".format(task.directory, file.name)
            # blob = storage.Blob(file.name, bucket)
            # Download the file
            with open(dest, 'wb') as file_obj:
                file.download_to_file(file_obj)

            # Untar the file
#            tar = tarfile.TarFile(dest, 'r:gz')
#            tar.extractall(os.path.dirname(dest))

    def _download_results_from_googlecloud(self):
        """Download all files from all instances/tasks to our result folder for femag

        """
        threads = []
        logger.info("Start downloading files")
        for t in self.job.tasks:
            thread = threading.Thread(target=self._download_bucket, args=(t, ))
            threads.append(thread)
            thread.start()

        self._wait_for_threads_finished(threads, "Downloading files")

    def _get_status_code(self, filename='exit_code'):
        """Get the status code from the caluclation
        Status code is written in a file

        Args:
            filename (str): The filename where the exit code is stored
        """
        status_code = []
        for t in self.job.tasks:
            dir = "{}/{}".format(t.directory, filename)
            file = open(dir, 'r')
            status_code.append(file.read())
        return status_code

    def _join(self, filename='exit_code', delete=True):
        """Wait until all calculation are finished:
        This means wait until every bucket has an file with the name of the given filename

        Args:
            filename (str): The filename where the exit code is stored
            delete (bool):  Should the instance be deleted after calculation
        """
        finished_tasks = []
        logger.info("Calculating..")
        gcs = storage.Client(self.project.project_id)
        while len(finished_tasks) < len(self.job.tasks):
            for t in [task for task in self.job.tasks if task not in finished_tasks]:
                bucket = gcs.get_bucket(t.id)
                if bucket.get_blob(filename):
                    finished_tasks.append(t)
                    logger.info("Calculation is finished for instance: {}".format(t.id))
                    if delete:
                        self._delete_instances([t])

            time.sleep(10)

        logger.info("Calculation is finished")

    def _cleanup(self):
        """Remove all buckets which belongs to this calculation round

        """
        threads = []
        for t in self.job.tasks:
            thread = threading.Thread(target=self._delete_bucket, args=(t.id, ))
            threads.append(thread)
            thread.start()

        logger.info("Deleting buckets: ")
        self._wait_for_threads_finished(threads, "Deleting buckets")

    def _delete_bucket(self, bucket_name):
        gcs = storage.Client(self.project.project_id)
        bucket = gcs.get_bucket(bucket_name)
        bucket.delete_blobs(bucket.list_blobs())
        bucket.delete(force=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # FEMAG STUFF
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def create_job(self, workdir):
        """Create a FEMAG :py:class:`CloudJob`

        Args:
            workdir (str): The workdir where the calculation files are stored

        Return:
            Cloud job (:py:class:`CloudJob`)
        """
        self.job = femagtools.job.CloudJob(workdir)
        return self.job

    def submit(self):
        """Starts the FEMAG calculation(s) on Google Cloud

        Return:
            length of started tasks (int)
        """
        self._create_data_buckets()
        self._upload_files_to_buckets()
        self._start_instances()
        return len(self.job.tasks)

    def join( self ):
        """Wait until all calculations are finished

        Return:
            list of all calculations status (C = Ok, X = error) (:obj:`list`)
        """
        self._join(filename='exit_code', delete=True)

        self._download_results_from_googlecloud()

        if int(self.config.get('DELETE_BUCKETS', 0)):
            self._cleanup()
        status = self._get_status_code("exit_code")
        for t, r in zip(self.job.tasks, status):
            t.status = 'C' if int(r)==0 else 'X'

        return status


if __name__ == "__main__":
    g_engine = Engine()
    pdb.set_trace()  # Set breakpoint
