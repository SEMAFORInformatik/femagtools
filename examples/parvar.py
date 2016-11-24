#!/usr/bin/env python
"""
 Parameter Variation with Femag
 """
import os
from femagtools.multiproc import Engine
# instead you can use on of the following
#
# from femagtools.condor import Engine
# from femagtools.amazon import Engine
# from femagtools.google import Engine
#
import femagtools.grid
import logging
import numpy as np

parvardef = {
    "objective_vars": [
        {"name": "dqPar.torque[-1]",
         "label": "Load Torque/Nm"},
        {"name": "torque[-1].ripple",
         "label": "Torque Ripple/Nm"},
        {"name": "machine.plfe[-1]",
         "label": "Iron Losses/W"}
    ],
    "population_size": 25,
    "decision_vars": [
        {"steps": 4, "bounds": [-50, 0],
         "name": "angl_i_up", "label":"Beta"},
        {"steps": 3, "bounds": [100, 200],
         "name": "current", "label":"Current/A"}
    ]
}

operatingConditions = {
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
    "calc_fe_loss": 1,
    "speed": 50.0,
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
    "relperm": 1.05,
    "magncond": 833333,
    "magnwidth": 15.0e-3,
    "magnlength": 100.0e-3,
    "hc_min": 760000.0}
]

magnetizingCurve = "./magnetcurves"

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
        "mcvkey_yoke": "M330-50A",
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
    "windings": {
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

engine = Engine()

userdir = os.path.expanduser('~')
workdir = os.path.join(userdir, 'parvar')
try:
    os.makedirs(workdir)
except OSError:
    pass

parvar = femagtools.grid.Grid(workdir,
                              magnetizingCurves=magnetizingCurve,
                              magnets=magnetMat)

results = parvar(parvardef, machine, operatingConditions, engine)

x = femagtools.grid.create_parameter_range(results['x'])
f = np.reshape(results['f'], (np.shape(results['f'])[0], np.shape(x)[0])).T

# print header
print(' '.join(['{:15}'.format(s)
                for s in [d['label']
                          for d in parvardef['decision_vars']] +
                [o['label']
                 for o in parvardef['objective_vars']]]))
print()
# print values in table format
for l in np.hstack((x, f)):
    print(' '.join(['{:15.2f}'.format(y) for y in l]))

# create scatter plot
#
import matplotlib.pyplot as pl
import mpl_toolkits.mplot3d as mpl
fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x[:, 0], x[:, 1], np.array(f[:, 0]))
ax.set_xlabel(parvardef['decision_vars'][0]['label'])
ax.set_ylabel(parvardef['decision_vars'][1]['label'])
ax.set_zlabel(parvardef['objective_vars'][0]['label'])
pl.savefig('parvar.png')
