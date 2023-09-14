import femagtools
import os
import logging


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
        num_slots_gen=3,
        mcvkey_yoke="dummy",
        rlength=1.0,
        stator1=dict(
            slot_rf1=0.057,
            tip_rh1=0.037,
            tip_rh2=0.037,
            tooth_width=0.009,
            slot_width=0.003)
    ),

    magnet=dict(
        mcvkey_shaft="dummy",
        mcvkey_yoke="dummy",
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

    winding=dict(
        num_phases=3,
        num_wires=100,
        coil_span=3.0,
        num_layers=1)
)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

workdir = os.path.join(
    os.path.expanduser('~'), 'femag')
try:
    os.makedirs(workdir)
except OSError:
    pass

femag = femagtools.Femag(workdir)

operatingConditions = dict(
    angl_i_up=[0.0, 0.0],
    calculationMode="pm_sym_loss",
    wind_temp=60.0,
    magn_temp=60.0,
    current=[50.0, 50.0],
    speed=[50.0, 100.0])

r = femag(machine,
          operatingConditions)

print("""
Speed   Torque  Current    Teeth    Yoke    Rotor    Total
1/min    Nm       A         W         W      W         W
""")
for x in zip(
        r['speed'],
        r['torque'],
        r['i1'],
        r['staza'],
        r['stajo'],
        r['rotfe'],
        r['total']):
    print("{0:6.1f}{1:9.1f}{2:8.2f}{3:9.1f}{4:9.1f}{5:8.1f}{6:10.1f}".format(
        60*x[0], x[1], x[2], x[3], x[4], x[5], x[6]))
