#!/usr/bin/env python
"""
 Parameter Variation with Femag
 """
import sys
import os
import json
import femagtools
import logging

opt = {
    "objective_vars": [
        {"desc": "Torque / Nm", "name": "dqPar.torque", "sign": -1},
        {"desc": "Torque Ripple / Nm", "name": "torque.ripple", "sign": 1},
        {"desc": "Iron Loss / W", "name": "machine.plfe", "sign": 1}
    ],
    "population_size": 3,
    "decision_vars": [
        {"desc": "beta", "steps": 3, "bounds": [-50, 0],
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

pmMotor = {
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
#engine = femagtools.MultiProc(numprocs)
engine = femagtools.Condor()

# Cloud computing
# engine = femagtools.DockerEngine(numprocs)

# precalculated AMAZONE S3
a_buckets = [{'id': '0422f2d4-0072-4eef-9ed9-01cb0a026c47.semafor.ch', 'folder': '/home/mau/parvar/0'},
             {'id': '1d859bf5-7a98-4218-9bda-ed72684e6b58.semafor.ch', 'folder': '/home/mau/parvar/1'},
             {'id': '820b5e48-91f0-458e-9e1a-d74447fc697a.semafor.ch', 'folder': '/home/mau/parvar/2'}]

# precalculated Google cloud storage
g_buckets = [{'id': '91d56074-c4c1-422a-99a2-1a0e5925b8e1-2', 'folder': '/home/mau/parvar/2'},
             {'id': 'fd29120d-686d-4213-9759-241c1c95dc5e-1', 'folder': '/home/mau/parvar/1'},
             {'id': '77e021b2-2f06-4b68-a9e8-5a862a6ff6ae-0', 'folder': '/home/mau/parvar/0'}]

# engine = femagtools.AmazonEngine(company='semafor',
#                                 number_of_instances=numprocs,
#                                 server_location='eu-central-1',
#                                 cluster_name='femag_semafor3',  # Cluster name
#                                 buckets=a_buckets,                # buckets with data
#                                 image_id='ami-b0cc23df',        # default amazon image for docker container
#                                 instance_type='t2.small',       # t2.small or t2.micro
#                                 task_definition='femag:10')     # femag:10 or femag:7

#engine = femagtools.GoogleEngine('semafor', buckets=g_buckets)

userdir = os.path.expanduser('~')
workdir = os.path.join(userdir, 'parvar')
try:
    os.makedirs(workdir)
except OSError:
    pass

grid = femagtools.Grid(workdir,
                       magnetizingCurves=magnetizingCurve,
                       magnets=magnetMat)

results = grid(opt, pmMotor, operatingConditions, engine)

json.dump(results, sys.stdout)
