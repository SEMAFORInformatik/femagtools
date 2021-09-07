#!/usr/bin/env python
"""
 Dakota Moat with Femag
 """
import pathlib
import json
from femagtools.multiproc import Engine
# instead you can use on of the following
#
#from femagtools.docker import Engine
#from femagtools.condor import Engine
# fr
# from femagtools.google import Engine
#
import femagtools.dakota

import logging
import numpy as np

parvardef = {
    "objective_vars": [
        {"name": "machine.torque",
         "label": "Load Torque/Nm"},
        {"name": "torque[0].ripple",
         "label": "Cogging Torque/Nm"},
        {"name": "torque[-1].ripple",
         "label": "Torque Ripple/Nm"}
    ],
    "population_size": 9,
    "decision_vars": [
        {"steps": 3, "bounds": [2e-3, 4e-3],
         "name": "stator.statorRotor3.slot_width",
         "label": "Slot Width/m"},
        {"steps": 3, "bounds": [0.75, 0.85],
         "name": "magnet.magnetSector.magn_width_pct",
         "label": "Rel. Magnet Width"},
        {"steps": 3, "bounds": [0.021, 0.0335],
         "name": "magnet.magnetSector.magn_shape",
         "label": "Magnet Shape/m"}
    ]
}

# Machine and Simulation Model

magnetMat = [{
    "name": "std",
    "remanenc": 1.2,
    "relperm": 1.05}]


machine = dict(
    name="PM-4-130",
    lfe=0.1,
    poles=4,
    outer_diam=0.13,
    bore_diam=0.07,
    inner_diam=0.015,
    airgap=0.0015,

    stator=dict(
        num_slots=12,
        num_slots_gen=3,
        mcvkey_yoke="M330-50A",
        rlength=1.0,
        statorRotor3=dict(
            slot_height=0.02,
            slot_h1=0.002,
            slot_h2=0.004,
            slot_r1=0.0,
            slot_r2=0.0,
            wedge_width1=0.0,
            wedge_width2=0.0,
            middle_line=0,
            tooth_width=0.009,
            slot_top_sh=0,
            slot_width=0.003)
    ),

    magnet=dict(
        material='std',
        mcvkey_yoke="M330-50A",
        magnetSector=dict(
            magn_num=1,
            magn_width_pct=0.6,
            magn_height=0.005,
            magn_shape=0.0,
            bridge_height=0,
            magn_type=2,
            condshaft_r=0.02,
            magn_ori=1,
            magn_rfe=0.0,
            bridge_width=0,
            magn_len=1)
    ),

    windings=dict(
        num_phases=3,
        num_wires=20,
        coil_span=3.0,
        num_layers=1)
)


simulation = dict(
    speed=5000.0 / 60,
    calculationMode="pm_sym_fast",
    magn_temp=20.0,
    wind_temp=60,
    period_frac=6,
    current=28.2842712474,
    angl_i_up=0.0)

magnetizingCurve = "../magnetcurves"


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

if __name__ == '__main__':
    workdir = pathlib.Path.home() / 'dakota'
    workdir.mkdir(parents=True, exist_ok=True)

    # try List, Grid, Sobol, LatinHypercube
    moat = femagtools.dakota.PsuadeMoat(workdir,
                                        magnetizingCurves=magnetizingCurve,
                                        magnets=magnetMat)

    # start calculation
    results = moat(parvardef, machine, simulation,
                   # engine=dict(module='femagtools.docker',
                   #            num_threads=8, port=5555
                   #            dispatcher='localhost'))
                   engine=dict(module='femagtools.multiproc',
                               process_count=8))

    with open('results.json', 'w') as fp:
        json.dump(results, fp)

    # print results
    x = results['x']
    f = results['f']

    # print header
    print(' '.join(['{:15}'.format(s)
                    for s in [d['label']
                              for d in parvardef['decision_vars']] +
                    [o['label']
                     for o in parvardef['objective_vars']]]))
    print()
    # print values in table format
    for l in np.vstack((x, f)).T:
        print(' '.join(['{:15.6f}'.format(y) for y in l]))
