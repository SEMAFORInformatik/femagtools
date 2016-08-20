#!/usr/bin/env python
"""
 Parameter Variation with Femag
 """
import os
import femagtools
import logging

parvar = {
    "objective_vars": [
        {"name": "dqPar.torque"},
        {"name": "torque.ripple"},
        {"name": "machine.plfe"}
    ],
    "population_size": 3,
    "decision_vars": [
        {"steps": 3, "bounds": [-50, 0],
         "name": "angl_i_up"}
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
        "zeroangle": 0.0,
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

numprocs = 3
engine = femagtools.MultiProc(numprocs)

userdir = os.path.expanduser('~')
workdir = os.path.join(userdir, 'parvar')
try:
    os.makedirs(workdir)
except OSError:
    pass

grid = femagtools.Grid(workdir,
                       magnetizingCurves=magnetizingCurve,
                       magnets=magnetMat)

results = grid(parvar, machine, operatingConditions, engine)

print("Beta    T/Nm  Tr/Nm   Pfe/W")
for l in zip(*(results['x']+results['f'])):
    print('{0:6.1f} {1:6.1f}{2:6.1f} {3:8.1f}'.format(*l))
