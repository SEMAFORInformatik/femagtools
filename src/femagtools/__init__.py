"""Python API for FEMAG

"""
__title__ = 'femagtools'
__version__ = '1.8.10'
__author__ = 'Ronald Tanner'
__license__ = 'BSD'
__copyright__ = 'Copyright 2023-2025 Gamma Technology'

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
    """create FSL command list from model and operating parameters

    Args:
        machine: dict with parameters (model)
        operatingConditions: dict with parameters (simulation)
        magnetmat: list fo dict with parameters (magnet material)
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
