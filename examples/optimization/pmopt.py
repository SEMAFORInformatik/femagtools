#!/usr/bin/env python
"""
 Multiobjective Optimization with Femag
 """
import sys
import json
import femagtools.opt
import logging
import glob
import pathlib
import os

from femagtools.multiproc import Engine
# instead you can use on of the following
#
#from femagtools.condor import Engine
# from femagtools.amazon import Engine
# from femagtools.google import Engine

opt = {
    "objective_vars": [
        {"desc": "Torque / Nm", "name": "dqPar.torque[-1]", "sign": -1},
        {"desc": "Torque Ripple / Nm", "name": "torque[-1].ripple"},
        {"desc": "Iron Loss / W", "name": "machine.plfe[-1]"}
    ],
    "population_size": 16,
    "decision_vars": [
        {"desc": "Magn width", "bounds": [0.75, 0.85],
         "name": "magnet.magnetSector.magn_width_pct"},
        {"desc": "Magn height", "bounds": [3e-3, 5e-3],
         "name": "magnet.magnetSector.magn_height"}
    ]
}

operatingConditions = {
    "angl_i_up": -10.0,
    "calculationMode": "pm_sym_fast",
    "wind_temp": 60.0,
    "magn_temp": 60.0,
    "current": 10.0,
    "eval_force": 0,
    "speed": 50.0,
    "period_frac": 6,
    "optim_i_up": 0
}

magnetMat = [{
    "name": "M395",
    "remanenc": 1.17,
    "temcoefbr": -0.001,
    "spmaweight": 7.5,
    "magntemp": 20.0,
    "temcoefhc": -0.001,
    "hcb": 810000.4,
    "relperm": 1.05}
]

magnetizingCurve = "../magnetcurves"

machine = dict(
    name="PM 130 L4",
    lfe=0.1,
    poles=4,
    outer_diam=0.13,
    bore_diam=0.07,
    inner_diam=0.015,
    airgap=0.001,

    stator=dict(
        num_slots=12,
        mcvkey_yoke="M270-35A",
        num_slots_gen=3,
        nodedist=1.5,
        rlength=1.0,
        statorRotor3=dict(
            slot_h1=0.002,
            slot_h2=0.004,
            middle_line=0,
            tooth_width=0.009,
            wedge_width2=0.0,
            wedge_width1=0.0,
            slot_top_sh=0,
            slot_r2=0.002,
            slot_height=0.02,
            slot_r1=0.003,
            slot_width=0.003)
    ),

    magnet=dict(
        nodedist=1.0,
        mcvkey_shaft="M270-35A",
        mcvkey_yoke="M270-35A",
        material="M395",
        magnetSector=dict(
            magn_num=1,
            magn_width_pct=0.8,
            magn_height=0.004,
            magn_shape=0.0,
            bridge_height=0.0,
            magn_type=1,
            condshaft_r=0.02,
            magn_ori=2,
            magn_rfe=0.0,
            bridge_width=0.0,
            magn_len=1.0)
    ),

    windings=dict(
        num_phases=3,
        num_wires=100,
        coil_span=3.0,
        slot_indul=1,
        cufilfact=0.45,
        culength=1.4,
        num_layers=1)
)


def get_losses(task):
    """adds extra results for optimization:
         total losses and efficiency"""
    bch = femagtools.bch.read(
        glob.glob(os.path.join(task.directory, '*.B*CH'))[-1])

    pltot = sum((bch.losses[-1]['winding'],
                 bch.losses[-1]['staza'], bch.losses[-1]['stajo'],
                 bch.losses[-1]['rotfe'],
                 bch.losses[-1]['magnetJ']))
    eff = bch.machine['p2']/(bch.machine['p2'] + pltot)

    bch.machine['pltotal'] = [pltot]
    bch.machine['eff'] = eff

    return bch


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    logger = logging.getLogger('pmopt')
    engine = Engine()

    workdir = pathlib.Path.home() / 'opti'
    workdir.mkdir(parents=True, exist_ok=True)

    o = femagtools.opt.Optimizer(workdir,
                                 magnetizingCurve, magnetMat,
                                 result_func=get_losses)

    num_generations = 3
    results = o.optimize(num_generations,
                         opt, machine, operatingConditions, engine)

    json.dump(results, sys.stdout)
