"""
   Parameter Identification of a IM
   Ronald Tanner
"""
import logging
import femagtools.machine
import femagtools.mcv
import pathlib
import json
from machine import machine, magnetizingCurves, condMat
from femagtools.multiproc import Engine

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    workdir = pathlib.Path('work')
    workdir.mkdir(exist_ok=True)

    engine = Engine()
    f1 = 50
    u1 = 400
    wdgcon = 'star'
    results = femagtools.machine.im.parident(
        workdir, engine, f1, u1, wdgcon,
        machine=machine,
        magnetizingCurves=magnetizingCurves,
        condMat=condMat)
    with open('impar.json', 'w') as fp:
        json.dump(results, fp)
