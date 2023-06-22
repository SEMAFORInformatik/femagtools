"""
  Parameter Identification of a EESM
  Ronald Tanner
"""
import logging
import json
import pathlib
from machine import machine, laminations, condMat
import femagtools.machine.sm
import femagtools.mcv
from femagtools.multiproc import Engine

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    workdir = pathlib.Path('work')
    workdir.mkdir(exist_ok=True)

    engine = Engine()
    r = femagtools.machine.sm.parident(workdir, engine, machine,
                                       magnetizingCurves=laminations,
                                       condMat=condMat)
    r['zeta1'] = 0.3  # skin effect factor of winding resistance
    with open('eecpars.json', 'w') as fp:
        json.dump(r, fp)
