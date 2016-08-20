# -*- coding: utf-8 -*-
"""
    femagtools
    ~~~~~~~~~~

    Python bindings for FEMAG

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
__title__ = 'femagtools'
__version__ = '0.0.1'
__author__ = 'Ronald Tanner'
__license__ = 'BSD'
__copyright__ = 'Copyright 2016 SEMAFOR Informatik & Energie AG'

from .bch import BchReader
from .machine import PmRelMachineLdq, PmRelMachinePsidq
from .femag import BaseFemag, Femag, ZmqFemag
from .model import Model, FeaModel
from .fsl import FslBuilder
from .magcurv import MagnetizingCurve
from .magnet import Magnet
from .grid import Grid
from .vbfreader import VbfReader
from .tksreader import TksReader
from .mcvreader import McvReader
from .poc import Poc
from .opt import Optimizer
from .condor import Condor
from .job import Job, CondorJob, AmazonJob
from .moproblem import FemagMoProblem
from .multiproc import MultiProc


def read_bchfile(filename):
    """Read BCH/BATCH results from file *filename*."""
    import io
    bchresults = BchReader()
    with io.open(filename) as f:
        bchresults.read(f.readlines())
    return bchresults


def create_fsl(machine, operatingconditions=None, magnetmat=None):
    """create FSL command list from model parameters"""
    model = Model(machine)
    builder = FslBuilder()
    if operatingconditions:
        if magnetmat:
            magnets = Magnet(magnetmat)
            return builder.create(model, operatingconditions, magnets)
        return builder.create(model, operatingconditions, None)
    return builder.create_model(model)
    
