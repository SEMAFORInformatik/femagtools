
Welcome to FEMAG Tools's documentation!
=======================================

Femagtools is an Open-Source Python-API for FEMAG offering following features:

* run Femag with a FSL file anywhere:
  locally (single and multi-core), remote (ZMQ), HT Condor, Cloud (Amazon AWS, Google Cloud)
* read BCH/BATCH files
* create FSL files from model and calculation templates
* calculate machine characteristics by using analytic machine models
* execute parameter studies and multi-objective optimization

Installation
------------

Install Femagtools with Pip::

  $ pip install femagtools

Prerequisite: a fairly recent FEMAG version (see http://www.profemag.ch) must be found in one of the
directories listed in your PATH variable.

Contents
========

.. toctree::
   :maxdepth: 2

   intro
   models
   bchreader
   amazon_engine
   google_engine
   apidocs


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

