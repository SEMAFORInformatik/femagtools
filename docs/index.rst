
Welcome to FEMAG Tools's documentation!
=======================================

   .. image:: img/femagtools.png


Femagtools is an Open-Source Python-API for FEMAG offering following features:

* run Femag with a FSL file anywhere:
  locally (single and multi-core), remote (ZMQ), HT Condor, Cloud (Amazon AWS, Google Cloud)
* read BCH/BATCH, I7/ISA7, PLT, ERG files
* read and write MCV files (magnetizing curves)
* create FSL files from model and calculation templates and/or user specific FSL or from DXF
* create a variety of plots
* calculate machine characteristics by using analytic machine models
* execute parameter studies and multi-objective optimization

The package is hosted on github: <https://github.com/SEMAFORInformatik/femagtools/>`.
Contributions are highly welcome.
  
Installation
------------

Femagtools can be installed on any 2.7 or 3.x Python distribution with Numpy, Scipy and Pip::

  $ pip install femagtools

Or with Anaconda::
  
  $ conda install -c semafor -c conda-forge femagtools

Prerequisite: a fairly recent FEMAG version
(see http://www.profemag.ch) must be found in one of the
directories listed in your PATH variable.

If a proxy is needed::

  $ pip --proxy http://proxy.hell:3128 install femagtools

For Anaconda you can setup the environment variables HTTP_PROXY and
HTTPS_PROXY or the file .condarc::

  proxy_servers:
    http: http://user:pass@corp.com:8080
    https: https://user:pass@corp.com:8080

For further information: <https://conda.io/docs/user-guide/configuration>_.

Console Scripts
===============

The following scripts can be executed from console:

* femagtools-plot: create plots from BCH/BATCH file
* femagtools-convert: various mesh format conversion
* femagtools-bchxml: convert BCH/BATCH file into XML
* femagtools-dxfsl: convert DXF into FSL


User Guide
==========

.. toctree::
   :maxdepth: 2

   intro
   femag
   models
   bchreader
   forcedens
   engine

References
==========

.. toctree::
   :maxdepth: 2

   apidocs


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Copyright
=========

Copyright: 2017-2018 Semafor Informatik & Energie AG, Basel

License: BSD, see LICENSE for more details.
