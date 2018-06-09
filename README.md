
# Introduction to Femagtools

![logo](https://github.com/SEMAFORInformatik/femagtools/raw/master/docs/img/femagtools.png)

[![Build Status](https://travis-ci.org/SEMAFORInformatik/femagtools.svg?branch=master)]
(https://travis-ci.org/SEMAFORInformatik/femagtools)

Femagtools is an Open-Source Python-API for FEMAG offering following features:

* run Femag with a FSL script file anywhere:
  locally (single and multi-core), remote (ZMQ), HT Condor, Cloud (Amazon AWS, Google Cloud)
* read I7/ISA7, BCH/BATCH, PLT files
* read and write MCV files (magnetizing curves)
* create a variety of plots
* create FSL files from model and calculation templates and/or user specific FSL
* create FSL files from DXF
* calculate machine characteristics by using analytic machine models
* execute parameter studies and multi-objective optimization

The package can be used with Python 2.7 or 3.x on Linux or Windows and is hosted on github: <https://github.com/SEMAFORInformatik/femagtools/>` where also many examples can be found in the examples directory. Contributions and feedback to this project are highly welcome.

The installation can be done in the usual ways with pip or conda.

For details see the documentation <http://docs.semafor.ch/femagtools>

## Modules and Scripts

The package provides following modules:

* __mcv__, __tks__, __jhb__, __losscoeffs__: handling magnetizing curves and iron losses
* __erg__, __bch__: read ERG, BCH/BATCH files created by FEMAG
* __model__, __fsl__: create machine and calculation models
* __femag__: manage the FEMAG calculation
* __airgap__: read airgap induction file created by a previous calculation
* __machine__: analytical machine models
* __grid__: running parameter variations
* __opt__: running multi objective optimizations
* __plot__: creating a variety of plots
* __dxfsl__: create FSL from DXF
* __isa7__: read ISA7/I7 files
* __forcedens__: read PLT files
* __amazon__, __google__, __condor__, __multiproc__: engines for the calculation in Cloud and HTCondor environments or locally using multiple cores

The following modules can be executed as script:

* __bch__: print content in json format if invoked witha BCH/BATCH file as argument
* __bchxml__: produces an XML file when invoked with a BCH/BATCH file as argument
* __plot__: produces a graphical report of a BCH/BATCH file
* __airgap__: prints the base harmonic amplitude of the radial component of the airgap induction when invoked with the file name of an airgap induction file
* __mcv__: print content in json format if invoked with a MC/MCV file as argument
* __dxfsl/conv__: show geometry or create fsl from dxf

## Usage
For many applications it is sufficient to import femagtools:


```python
import femagtools
```

The version can be checked with:


```python
femagtools.__version__
```




    '0.4'


