#!/usr/bin/env python
"""
 Parameter Variation with Femag
 """
import pathlib
import json
#from femagtools.multiproc import Engine
# instead you can use on of the following
#
from femagtools.docker import Engine
#from femagtools.condor import Engine
# fr
# from femagtools.google import Engine
#
# chose sampling method
import femagtools.parstudy

import femagtools.poc
import logging
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D

parvardef = {
    "objective_vars": [
        {"name": "machine.torque",
         "label": "Load Torque/Nm"},
        {"name": "torque[-1].ripple",
         "label": "Torque Ripple/Nm"},
        {"name": "machine.plfe[-1]",
         "label": "Iron Losses/W"},
        {"name": "airgap.Bamp",
         "label": "Airgap Induction/T"}
    ],
    "population_size": 25,
    "decision_vars": [
        {"steps": 3, "bounds": [-50, 0],
         "name": "angl_i_up", "label":"Beta"},
        {"steps": 3, "bounds": [100, 200],
         "name": "current", "label":"Current/A"}
    ]
}

poc = femagtools.poc.HspPoc(harm=[1, 5],
                            amp=[1, 0.01],
                            phi=[0, 0])

simulation = {
    "num_move_steps": 49,
    "angl_i_up": 0.0,
    "calculationMode": "pm_sym_fast",
    "wind_temp": 60.0,
    "magn_temp": 60.0,
    "current": 250.0,
    "eval_force": 0,
    "skew_angle": 0.0,
    "num_par_wdgs": 1,
    "num_skew_steps": 0,
    "period_frac": 6,
    "calc_fe_loss": 1,
    "speed": 50.0,
    "poc": poc,
    "optim_i_up": 0,
    "airgap_induc": True
}

magnetMat = [{
    "name": "M395",
    "remanenc": 1.17,
    "temcoefbr": -0.001,
    "spmaweight": 7.5,
    "magntemp": 20.0,
    "temcoefhc": -0.001,
    "HcB": 810000.4,   # Note: relperm = remanenc/mue0/HcB
    "relperm": 1.05,
    "magncond": 833333,
    "magnwidth": 15.0e-3,
    "magnlength": 100.0e-3,
    "HcJ": 760000.0}
]

magnetizingCurve = "../magnetcurves"

machine = {
    "name": "PM 270 L8",
    "desc": "PM Motor 270mm 8 poles VMAGN",
    "poles": 8,
    "outer_diam": 0.26924,
    "bore_diam": 0.16192,
    "inner_diam": 0.11064,
    "airgap": 0.00075,
    "lfe": 0.08356,
    "stator": {
        "num_slots": 48,
        "num_slots_gen": 12,
        "fillfac": 0.96,
        "mcvkey_yoke": "M330-50A",
        "mcvkey_teeth": "M270-35A",
        "nodedist": 4.0,
        "statorRotor3": {
            "slot_height": 0.0335,
            "slot_h1": 0.001,
            "slot_h2": 0.0,
            "slot_width": 0.00193,
            "slot_r1": 0.0001,
            "slot_r2": 0.00282,
            "wedge_width1": 0.00295,
            "wedge_width2": 0.0,
            "middle_line": 0.0,
            "tooth_width": 0.0,
            "slot_top_sh": 0.0}
    },
    "magnet": {
        "nodedist": 1.0,
        "material": "M395",
        "mcvkey_yoke": "M330-50A",
        "magnetIronV": {
            "magn_angle": 145.0,
            "magn_height": 0.00648,
            "magn_width": 0.018,
            "condshaft_r": 0.05532,
            "magn_num": 1.0,
            "air_triangle": 1,
            "iron_hs": 0.0001,
            "gap_ma_iron": 0.0002,
            "iron_height": 0.00261,
            "magn_rem": 1.2,
            "iron_shape": 0.0802
        }
    },
    "winding": {
        "num_phases": 3,
        "num_layers": 1,
        "num_wires": 9,
        "coil_span": 6.0,
        "cufilfact": 0.4,
        "culength": 1.4,
        "slot_indul": 0.5e-3
    }
}

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

if __name__ == '__main__':
    engine = Engine(num_threads=8)

    workdir = pathlib.Path.home() / 'parvar'
    workdir.mkdir(parents=True, exist_ok=True)

    # try List, Grid, Sobol, LatinHypercube
    parvar = femagtools.parstudy.LatinHypercube(workdir,
                                                magnetizingCurves=magnetizingCurve,
                                                magnets=magnetMat)

    # keep the BCH/BATCH files of each run
    repdir = workdir / 'report'
    workdir.mkdir(parents=True, exist_ok=True)
    # parvar.set_report_directory(repdir)

    # start calculation
    results = parvar(parvardef, machine, simulation,
                     engine, num_samples=8)

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
        print(' '.join(['{:15.2f}'.format(y) for y in l]))

    # create scatter plot
    #
    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x[0], x[1], np.array(f[0]))
    ax.set_xlabel(parvardef['decision_vars'][0]['label'])
    ax.set_ylabel(parvardef['decision_vars'][1]['label'])
    ax.set_zlabel(parvardef['objective_vars'][0]['label'])
    pl.savefig('parvar.png')
