# Axial Flux Machine (AFM) demo
#
from machine import machine

if __name__ == '__main__':
    import logging
    import pathlib
    import json
    from femagtools.machine.afpm import parident
    from femagtools.multiproc import Engine

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    workdir = pathlib.Path('dqpar')
    workdir.mkdir(exist_ok=True)

    engine = Engine()
    temp = [60, 90]
    r = parident(workdir, engine, temp, machine,
                 magnetizingCurves='',
                 magnetMat=[], condMat=[],
                 speed=4000/60,
                 beta_min=-90,
                 num_beta_steps=4, num_cur_steps=3)
    with open('dqpar.json', 'w') as fp:
        json.dump(r, fp)
