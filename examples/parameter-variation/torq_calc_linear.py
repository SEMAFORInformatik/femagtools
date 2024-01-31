import pathlib
import logging
import json
import numpy as np
import femagtools.parstudy
from femagtools.multiproc import Engine

poc = femagtools.poc.Poc(50, parameters=dict(
    phi_voltage_winding=[251, 11, 131]))
# D-Axis: phi_voltage_winding[0] - 90 --> 161

machine = dict(
    name="PM lin",
    lfe=0.05,
    poles=10,
    airgap=0.0015,
    coord_system=1,

    stator=dict(
        num_slots=12,
        num_slots_gen=6,
        mcvkey_yoke="dummy",
        rlength=1.0,
        stator3Linear=dict(
            slot_height=0.02,
            slot_h1=0.002,
            slot_h2=0.002,
            tip_slot=0.003,
            yoke_height=0.008,
            slot_r1=0.004,
            slot_r2=0.005,
            tooth_width=0.01,
            width_bz=0.025,
            middle_line=1)
    ),
    magnet=dict(
        mcvkey_yoke="dummy",
        magnetSectorLinear=dict(
            magn_height=0.008,
            magn_width_pct=0.8,
            pole_width=0.03,  # bz * Q/P
            yoke_height=0.008,
            magn_len=1.0,
            gap_ma_yoke=0,
            magn_ori=0,
            airgap_shape=0.0,
            magn_type=1)
    ),
    winding=dict(
        num_phases=3,
        num_wires=20,
        coil_span=1.0,
        num_layers=2)
)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

parvardef = {
    "objective_vars": [
        {"name": "linearForce_fft[0].force[0]",
         "label": "Force / N"}],
    "decision_vars": [
        {"steps": 3, "bounds": [-50, 0],
         "name": "angl_i_up", "label":"Beta"},
        {"steps": 3, "bounds": [0, 10],
         "name": "current", "label":"Current/A"}
    ]
}

simulation = dict(
    angl_i_up=0.0,
    calculationMode="torq_calc",
    wind_temp=20.0,
    magn_temp=20.0,
    current=7.07,
    poc=poc,
    speed=10.0)


engine = Engine(num_threads=8)

workdir = pathlib.Path.home() / 'parvar'
workdir.mkdir(parents=True, exist_ok=True)

# try List, Grid, Sobol, LatinHypercube
parvar = femagtools.parstudy.LatinHypercube(workdir)
#                                            magnetizingCurves=magnetizingCurve,
#                                            magnets=magnetMat)

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
