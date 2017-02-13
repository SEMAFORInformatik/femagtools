#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Ld-Lq-Identification with Femag
 """
import os
import femagtools
import femagtools.machine
import logging
import numpy as np

feapars = {
    "num_move_steps": 25,
    "calculationMode": "psd_psq_fast",
    "magn_temp": 60.0,
    "maxid": 0.0,
    "minid": -120.0,
    "maxiq": 120.0,
    "miniq": -0.0,
    "delta_id": 40.0,
    "delta_iq": 40.0,
    "speed": 50.0
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

workdir = os.path.join(os.path.expanduser('~'), 'femag')
try:
    os.makedirs(workdir)
except OSError:
    pass

femag = femagtools.Femag(workdir,
                         magnetizingCurves=magnetizingCurve,
                         magnets=magnetMat)

r = femag(pmMotor, feapars)

print(r.type)

# find speed at u1max
u1max = 340
tq = 170

psid = r.psidq['psid']
psiq = r.psidq['psiq']
id = r.psidq['id']
iq = r.psidq['iq']

p = r.machine['p']
r1 = 0.0

pm = femagtools.machine.PmRelMachinePsidq(3, p,
                                          psid,
                                          psiq,
                                          r1,
                                          id,
                                          iq)

tq = 170.0
u1 = 340.0

iqx, idx = pm.iqd_torque(tq)
w1 = pm.w1_u(u1, idx, iqx)
i1 = np.linalg.norm(np.array((iqx, idx)))
betaopt = np.arctan2(idx, iqx)/np.pi*180

print("f1 {0:8.1f} Hz,  I1 {1:8.1f} A, Beta {2:4.1f} °".format(
    w1/2/np.pi, i1, betaopt))
