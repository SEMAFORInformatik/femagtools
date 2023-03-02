# -*- coding: utf-8 -*-
"""
    femagtools
    ~~~~~~~~~~

    Python bindings for FEMAG



"""
__title__ = 'femagtools'
__version__ = '1.2.4'
__author__ = 'Ronald Tanner'
__license__ = 'BSD'
__copyright__ = 'Copyright 2016-2022 SEMAFOR Informatik & Energie AG'

from .asm import read
from .bch import Reader
from .model import MachineModel
from .fsl import Builder
from .magnet import Magnet
from .femag import Femag, ZmqFemag
from .conductor import Conductor


def read_bchfile(filename):
    """Read BCH/BATCH results from file *filename*."""
    import io
    bchresults = Reader()
    with io.open(filename, encoding='latin1', errors='ignore') as f:
        bchresults.read(f.readlines())
    return bchresults


def create_fsl(machine,
               operatingconditions={},
               magnetmat=[],
               condMat=[],
               templatedirs=[]):
    """create FSL command list from model parameters

    Args:
        machine: dict with parameters
        operatuingConditions: dict with parameters
        magnetmat: list fo dict with parameters
        templatedir: additional template directory
"""
    model = MachineModel(machine)

    builder = Builder(templatedirs)
    magnets = []
    if magnetmat:
        magnets = Magnet(magnetmat)
    if condMat:
        if not isinstance(condMat, Conductor):
            condMat = Conductor(condMat)

    if operatingconditions:
        return builder.create(model, operatingconditions, magnets)
    return builder.create_model(model,
                                condMat=condMat,
                                magnets=magnets) + ['save_model("cont")']
