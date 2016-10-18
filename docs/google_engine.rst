Google Engine
*************

The Google engine calculates the FEMAG tasks in the Google Cloud.
The files for calculation are uploaded to an unique Google bucket.

Prerequisite
============
To use Google as a calculation engine you have to:

* create an account and setup a user with the correct authority.
* configured your google credentials with the google command tool

Fast start
==========
To use Google Cloud as engine, setup a configuartion and a startup bash file and use google as engine::

  import femagtools
  engine = femagtools.google.Engine()
  # Load config file
  engine.config.from_ini_file('config.ini')

If you already uploaded some FEMAG files you can add the buckets as list::

  buckets = [{'id': '926c39a7-db69-42f8-bf86-2fd1f8559188-0', 'folder': '~/parvar/0'},
             {'id': '99ec9af5-ecb5-4029-88d6-70da9841ef91-1', 'folder': '~/parvar/1'},
             {'id': '4c61ebb3-89df-4e3a-ba46-288f8e0672a6-2', 'folder': '~/parvar/2'}]
  
Configuration file
------------------
You have to setup a configuration file with an google Section, where you define some Options for the Google instance.

====================  =========================  =======================================================
Name                  Example value              Description
====================  =========================  =======================================================
server_location       us-east1-b                 Server location of your instances
instance_type         n1-standard-1              Google instance type
image_id*             femag                      The image of the instance. It is very advisable to use a custom image where femag and his dependencies are installed.
startup               ./startup.sh               A path to a startup.sh file which will be run when the instance is started. Usually here you define your calculation for femag
delete_buckets        1                          Delete the buckets after calculation is finished. 1 Delete 0 or no entry means no delete
server_config         ./gcloud.json              A json file where you can define some attributes for the instances
====================  =========================  =======================================================

All configuration options are optional.

Example of a configuration file::
  
 [google]
 SERVER_LOCATION = us-east1-b
 instance_type = n1-standard-1
 image_id =  femag-v1
 startup = ./cloud_init.txt
 company_name = Semafor
 delete_buckets = 1
 server_config = ./gcloud.json
 startup_script = ./startup.sh

 
startup.sh
----------
With the startup.sh file you can define, what happens, after the instance startet up. This is a good place to start the calculation.

.. note:: The files for one calculation are transferd in a tar.gz file

Example::
 
 #!/bin/bash

 # [START startup_script]
 CS_BUCKET=$(curl http://metadata/computeMetadata/v1/instance/attributes/bucket -H "Metadata-Flavor: Google")
 
 mkdir data
 gsutil rsync gs://$CS_BUCKET data
 cd data
 tar -xzf *.tar.gz
 xfemag -b femag.fsl
 echo $? > exit_code
 gsutil rsync /data gs://$CS_BUCKET
 
 # [END startup_script]

Server configuration
--------------------
The Google Cloud allows you to set the configuartion over a JSON file. It is recommendable to set only attributes which are valid for all instances.

These attributes are automaticly set during runtime if they are missing in the json configuration:

=============  ================
Value          Descritpion
=============  ================
sourceImage    Your project id with your image name
machineType    Your server location with your instance type
name           Same name as the buckets with your company name first
bucket_id      The bucket id to download the correct data bucket
=============  ================

If you want more data in your instance, please set them as metadata items, and load them in the instance.

Example::

  {
    "disks": [
        {
            "boot": true,
            "autoDelete": true
        }
    ],
    "networkInterfaces": [{
        "network": "global/networks/default",
        "accessConfigs": [
            {"type": "ONE_TO_ONE_NAT", "name": "External NAT"}
        ]
    }],
    "serviceAccounts": [{
        "email": "default",
        "scopes": [
            "https://www.googleapis.com/auth/devstorage.read_write",
            "https://www.googleapis.com/auth/logging.write"
        ]
    }],


    "metadata": {
        "items": [{
            "key": "text",
            "value": "Run femag"
        }]
    }
  }

