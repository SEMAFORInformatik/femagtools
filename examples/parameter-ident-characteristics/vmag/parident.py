"""
  DQ Parameter Identification of a VMAG PM
  Ronald Tanner
"""
import logging
import pathlib
import json
import femagtools.machine
import femagtools.mcv
from machine import machine, laminations, magnetMat, condMat
from femagtools.multiproc import Engine

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    workdir = pathlib.Path('work')
    workdir.mkdir(exist_ok=True)

    engine = Engine()
    dqpars = femagtools.machine.pm.parident(
        workdir, engine, temp=[60, 90],
        i1_max=300,
        machine=machine,
        magnetizingCurves=laminations,
        magnetMat=magnetMat, condMat=condMat)
    dqpars['zeta1'] = 0.3  # skin effect factor of winding resistance
    with open('dqpar.json', 'w') as fp:
        json.dump(dqpars, fp)
