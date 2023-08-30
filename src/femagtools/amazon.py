"""
    femagtools.amazon
    ~~~~~~~~~~~~~~~~~

    Running FEMAG on Amazon Cloud EC2




    .. note: To use this engine you have to install the boto3 module from amazon
"""
import os
import threading
import time
import logging

import femagtools.job
from .config import Config

logger = logging.getLogger(__name__)

class MissingConfigurationException(Exception):
    def __init__(self, message):
        Exception.__init__(self, "Missing configuration: {}".format(message))


class Engine(object):

    config_class = Config

    default_config = {
        'ENGINE': 'amazon',
        'SERVER_LOCATION': 'eu-central-1',
        'INSTANCE_TYPE': 't2.micro',
        'ACL': 'authenticated-read',
        'IMAGE_ID': 'ami-b0cc23df',
        'FINISH_TASK_FILENAME': 'exit_code',
        'COMPANY_NAME': 'femag'
    }

    """The Amazon Engine

    This engine uses the boto3 Python module to interact
       with the amazon ec2 and s3 services

    Args:
        buckets (:obj:`list`): Existing buckets with femag calculation files
        configfile (str): Filename of config file

    .. :note: If possible you should use the same location for all services

    """
    def __init__(self, buckets=None, configfile='config.ini'):
        self.buckets = buckets
        self.job = None

        # Amazon file storage
        self.s3_resource = self._create_amazon_resource('s3')
        # Amazon Server administration
        self.ec2_resource = self._create_amazon_resource('ec2')

        # Create instance of config
        self.config = Config(self.default_config)
        self.config.from_ini_file(configfile)

    def _create_amazon_resource(self, resource):
        import boto3
        return boto3.resource(resource)

    def _create_data_buckets(self):
        """Create unique S3 Buckets for calculation

        Args:
            ACL (str): ACL-Rules for Amazon
        """
        # If buckets exsists map them with a folder
        if self.buckets:
            for idx, bucket in enumerate(self.buckets):
                self.job.tasks[idx].id = bucket['id']
                self.job.tasks[idx].directory = bucket['folder']
            return

        bucketConfiguration = {'LocationConstraint': self.config['SERVER_LOCATION']}

        # Create a bucket for every calculation
        for t in self.job.tasks:
            self.s3_resource.create_bucket(ACL=self.config['ACL'],
                                           Bucket=t.id,
                                           CreateBucketConfiguration=bucketConfiguration)

        logger.debug("Created buckets")

    def _upload_files_to_s3(self):
        """Upload all files to Amazon S3 for this calculation

        """
        if self.buckets:
            logger.info("Files are already uploaded")
            return

        threads = []
        for t in self.job.tasks:
            thread = threading.Thread(target=self._upload, args=(t, ))
            threads.append(thread)
            thread.start()

        logger.info("Uploading files: ")
        self._wait_for_threads_finished(threads, "Upload files")

    def _upload(self, task):
        """Upload thread for uploading one directory

        :internal:

        Args:
            task (py:class:`CloudTask`): The task which belongs to the uploading folder
        """
        # Upload one single tar_file
        task.tar_file.close()

        name = os.path.basename(task.file)
        Body = open(task.file, 'rb')
        self.s3_resource.Object(task.id, name).put(Body=Body)

    def _wait_for_threads_finished(self, threads, operation):
        """Wait until all threads are finished

        :internal:

        Args:
            threads (:obj:`list`): List of threads to check if they are finished
            operation (str): Name of the operation to write a meaningful log message

        """
        # Wait until all threads are not running
        while not all([t.isAlive() is False for t in threads]):
            time.sleep(5)

        # timer.cancel
        logger.info("{} is finished".format(operation))

    def _start_instances(self):
        """Start all instances for the calculation

        """
        # Prepare arguemtns for instance start
        param = {'MinCount': 1, 'MaxCount': 1 }

        if self.config.get('IMAGE_ID', None):
            param['ImageId'] = self.config['IMAGE_ID']
        else:
            raise MissingConfigurationException('image_id')

        if self.config.get('INSTANCE_TYPE', None):
            param['InstanceType'] = self.config['INSTANCE_TYPE']
        else:
            raise MissingConfigurationException('instance_type')

        if self.config.get('IAM_INSTANCE_PROFILE', None):
            param['IamInstanceProfile'] = {'Name': self.config['IAM_INSTANCE_PROFILE'] }

        if self.config.get('KEY_NAME', None):
            param['KeyName'] = self.config['KEY_NAME']

        # Set security group id as list
        if self.config['SECURITY_GROUP_IDS']:
            param['SecurityGroupIds'] = []
            for security_group in [s for s in self.config['SECURITY_GROUP_IDS'] if s]:
                param['SecurityGroupIds'].append(security_group)

        param['DryRun'] = self.config.get('DRY_RUN', False)

        threads = []
        for idx, t in enumerate(self.job.tasks):
            thread = threading.Thread(target=self._start_instance, args=(param, t))
            threads.append(thread)
            thread.start()

        self._wait_for_threads_finished(threads, "Start instances")

    def _start_instance(self, param, task):
        """Start one instance

        :internal:

        Args:
            task (Task): the task for calculation

        """
        user_data = self._read_cloud_init(task.id)
        if user_data:
            param['UserData'] = user_data

        instances = self.ec2_resource.create_instances(**param)

        instance = instances[0]  # We only started one instance
        logger.info("Instance started: {}".format(instance.id))
        self._add_tag(task.id, instance.id)

        instance.wait_until_running()
        instance.load()  # Reload the data to get public dns etc.
        logger.info("Instance {} is running: Public dns: {}".format(instance.id, instance.public_dns_name))
        task.ec2_instance = instance.id

    def _add_tag(self, task_id, instance_id):
        """Add a tag to the instance

        :internal:

        Args:
            task_id (int): The task id (Same as the S3 Bucket name)
            instance_id (int): The instance_id to set the tag to the right instance
        """
        tag = '{}-{}'.format(task_id, self.config.get('COMPANY_NAME', 'femag'))
        self.ec2_resource.create_tags(Resources=[instance_id], Tags=[{'Key': 'Name', 'Value': tag}])

    def _read_cloud_init(self, bucket_name):
        """Read the cloud init file and if there is a line which starts with {{ENV}}
        then put all config options as environment variables.

        """
        user_data = ""
        # Set all config options as environment variable
        if os.path.isfile(self.config.get('CLOUD_INIT', None)):
            with open(self.config['CLOUD_INIT'], 'rt') as f:
                for line in f:
                    if line.startswith('{{ENV}}'):
                        # Add config
                        for key, value in sorted(self.config.items()):
                            user_data += "export {}={}\n".format(key, value)
                        # add other important stuff
                        user_data += "export BUCKET_NAME={}\n".format(bucket_name)
                        continue
                    user_data += line
        return user_data

    def _join(self, timeout=20, filename='exit_code'):
        """Wait until all instances are finished with the calulation.

        :internal:

        Args:
            timeout (int): How long we wait between a check
            filename (str): What is the filename of the exit_code
        """
        import botocore       # For exception

        finished_tasks = []
        client = self.s3_resource.meta.client
        while len(finished_tasks) < len(self.job.tasks):
            for t in [task for task in self.job.tasks if task not in finished_tasks]:
                try:
                    client.get_object(Bucket=t.id, Key=filename)
                except botocore.exceptions.ClientError:
                    # Instance not ready
                    time.sleep(2)
                    continue
                finished_tasks.append(t)
                logger.info("Calculation is finished for instance {}".format(t.id))
                self.ec2_resource.instances.filter(InstanceIds=[t.ec2_instance]).terminate()
            time.sleep(timeout)

        logger.info("Calculations are finished")

    def _get_result_data_from_S3(self):
        """Get all the calculated files to the correct folder

        """
        import boto3
        client = self.s3_resource.meta.client
        transfer = boto3.s3.transfer.S3Transfer(client)
        for t in self.job.tasks:
            bucket = t.id
            folder = t.directory
            files = client.list_objects(Bucket=bucket)['Contents']
            logger.debug("Starting new folder")

            for file in files:
                file_name = file['Key']
                transfer.download_file(bucket, file_name, os.path.join("{}/{}".format(folder, file_name)))
                logger.debug("Downloaded file {}".format(file_name))

    def _get_status_code(self, filename='exit_code'):
        """Get the status code from the caluclation

        Args:
            filename (str): Filename of exit_code
        """
        status_code = []
        for t in self.job.tasks:
            dir = "{}/{}".format(t.directory, filename)
            file = open(dir, 'r')
            status_code.append(file.read())
        return status_code

    def _cleanup(self):
        threads = []
        for t in self.job.tasks:
            thread = threading.Thread(target=self._delete_bucket, args=(t.id, ))
            threads.append(thread)
            thread.start()

        logger.info("Deleting buckets: ")
        self._wait_for_threads_finished(threads, "Deleting buckets")

        # Clean up volumes
        client = self.ec2_resource.meta.client
        volumes = client.describe_volumes(Filters=[{'Name': 'status', 'Values': ['available']}])['Volumes']
        for v in volumes:
            client.delete_volume(VolumeId=v['VolumeId'])

    def _delete_bucket(self, bucket_name):
        bucket = self.s3_resource.Bucket(bucket_name)
        for key in bucket.objects.all():
            key.delete()
        bucket.delete()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # FEMAG STUFF
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def create_job(self, workdir):
        """Create a FEMAG :py:class:`CloudJob`

        Args:
            workdir (str): The workdir where the calculation files are stored

        Return:
            Cloud job (:class:`CloudJob`)
        """
        self.job = femagtools.job.CloudJob(workdir)
        return self.job

    def submit(self):
        """Starts the FEMAG calculation(s) on Amazon

        Return:
            length of started tasks (int)
        """
        self._create_data_buckets()
        self._upload_files_to_s3()
        self._start_instances()
        return len(self.job.tasks)

    def join( self ):
        """Wait until all calculations are finished

        Return:
            list of all calculations status (C = Ok, X = error) (:obj:`list`)
        """
        status = []
        # Wait until all tasks are finished
        self._join(timeout=20, filename=self.config['FINISH_TASK_FILENAME'])

        # get all files
        self._get_result_data_from_S3()

        # Remove buckets if cleanup is set
        if int(self.config.get('DELETE_BUCKETS', 0)):
            self._cleanup()

        status = self._get_status_code(filename=self.config['FINISH_TASK_FILENAME'])
        for t, r in zip(self.job.tasks, status):
            t.status = 'C' if int(r)==0 else 'X'

        return status
