
Welcome to FEMAG Tools's documentation!
=======================================

   .. image:: img/femagtools.png


Femagtools is an Open-Source Python-API for FEMAG offering following features:

* run Femag with a FSL file anywhere:
  locally (single and multi-core), remote (ZMQ), HT Condor, Cloud (Amazon AWS, Google Cloud)
* read BCH/BATCH, I7/ISA7, NC, PLT, ERG files
* read and write MCV files (magnetizing curves)
* create FSL files from model and calculation templates and/or user specific FSL or from DXF
* create and analyse symmetrical windings
* create a machine and winding layout from basic requirements
* identify the model parameters of equivalent circuit
* calculate machine characteristics and efficiency maps by using analytic machine models
* create a variety of plots
* execute parameter studies and multi-objective optimization

A couple of jupyter notebooks and example scripts can be found in the source directory
hosted on github: <https://github.com/SEMAFORInformatik/femagtools/>`.

Contributions and feedback are highly welcome.

Installation
------------

Femagtools can be installed on any 3.x Python distribution with Numpy, Scipy and Pip::

  $ pip install femagtools

Note::

  The package has a couple of optional dependencies:

  * dxfsl = ['dxfgrabber', 'networkx']
  * mplot =  ['matplotlib']
  * meshio = ['meshio']
  * vtk = ['vtk']
  * zmq = ['pyzmq']

  If you want them all:

    $ pip install femagtools[all]

  or individually by appending the name of the dependency in brackets. Example:

    $ pip install femagtools[dxfsl]

Prerequisite: a fairly recent FEMAG version
(see <https://www.gtisoft.com/download/femag-download>) must be found in one of the
directories listed in your PATH variable.

If a proxy is needed::

  $ pip --proxy http://proxy.hell:3128 install femagtools

For further information: <https://pip.pypa.io>.

Console Scripts
---------------

The following scripts can be executed from console:

* femagtools-plot: create plots from BCH/BATCH file
* femagtools-convert: various mesh format conversion
* femagtools-bchxml: convert BCH/BATCH file into XML
* femagtools-dxfsl: convert DXF into FSL


User Guide
----------

.. toctree::
   :maxdepth: 2

   intro
   femag
   models
   sizing
   windings
   bchreader
   ncisa
   forcedens
   dakota
   engine
   amela
   tspost
..
   API reference <_autosummary/femagtools>

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Copyright
---------

Copyright:

* 2017-2023 Semafor Informatik & Energie AG, Basel
* 2024 Gamma Technology LLC, Westmont, Illinois

License: BSD, see LICENSE for more details.
