Amazon Engine
*************

The Amazon engine calculates the FEMAG tasks in the Amazon EC2 cloud.
The files for calculation are uploaded to an unique S3 Bucket.

Prerequisite
============
To use Amazon as a calculation engine you have to:

* create an account and setup a user with the correct authority.
* configured your aws credentials under ~/.aws/credentials

Fast start
==========
To use Amazon as engine, setup a configuration and a cloud init file and use amazon as engine::

  import femagtools
  engine = femagtools.amazon.Engine(buckets=None, configfile='config.ini')

If you already uploaded some FEMAG files you can add the bucket as list::

  buckets = [{'id': '926c39a7-db69-42f8-bf86-2fd1f8559188-0', 'folder': '~/parvar/0'},
             {'id': '99ec9af5-ecb5-4029-88d6-70da9841ef91-1', 'folder': '~/parvar/1'},
             {'id': '4c61ebb3-89df-4e3a-ba46-288f8e0672a6-2', 'folder': '~/parvar/2'}]
  
Configuration file
------------------
You have to setup a configuration file with an amazon Section, where you define some Amazon EC2 Instance attribute.

====================  =========================  =======================================================
Name                  Example value              Description
====================  =========================  =======================================================
iam_instance_profile  ecsInstanceRole            An Amazon Role name. The Role should at least have the AmazonEC2ContainerServiceforEC2Role and AmazonS3FullAccess to create EC2 Instances and S3 Buckets
key_name              test-key-pair-eucentral-1  SSH Key name for login over SSH. It is not neccessary but useful if something does not work as expected on the instance.
security_group_ids    sg-8475b5dd,               A list of security group ids. If you only have one security group it is important to add a , at the end.
server_location       eu-central-1               Server location of your instances
instance_type         t2.small                   Instance type. Femag does not run if you use an instance type under the t2.small
image_id*              ami-6d4ceb03               The image of the instance. It is very advisable to use a custom image where femag and his dependencies are installed.
cloud_init            ./cloud_init.txt           A path to a cloud_init file which will be run when the instance is started. Usually here you define your calculation for femag
delete_buckets        1                          Delete the buckets after calculation is finished. 1 Delete 0 or no entry means no delete
====================  =========================  =======================================================

All configuration options are optional.

Example of a configuration file::
  
  [amazon]
  iam_instance_profile = ecsInstanceRole
  key_name = test-key-pair-eucentral-1
  security_group_ids = sg-8475b5dd,
  server_location = eu-central-1
  instance_type = t2.small
  image_id =  ami-6d4ceb03
  cloud_init = ./cloud_init.txt

.. note:: The **security_group_ids** is a list. If you only have one security id you have to put a ',' at the end

Cloud Init
----------
With the cloud init file you can define the command to be excuted after the instance has started up. This is a good place to start the calculation.

The special entry *{{ENV}}* indicates the femagtools module to put all the configuration options as environment variables. Additionally it adds the bucket name as *BUCKET_NAME* environment to download the files from the correct Bucket.

.. note:: The files are for one calculation transferd as in a tar.gz file
Example::
 
 #!/bin/bash
 
 {{ENV}}
 
 mkdir ~/.aws
 echo -e '[default]\nregion = '$SERVER_LOCATION'\noutput = json' > ~/.aws/config
 
 yum install -y aws-cli libquadmath
 
 mkdir ~/data
 aws s3 sync s3://$BUCKET_NAME/ ~/data
 cd ~/data
 tar -xzf *.tar.gz
 /usr/local/bin/xfemag -b femag.fsl </dev/null
 echo $? > ~/data/exit_code
 aws s3 sync ~/data s3://$BUCKET_NAME
