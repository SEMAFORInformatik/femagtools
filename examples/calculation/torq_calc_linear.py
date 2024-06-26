import femagtools
import os
import logging

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

workdir = os.path.join(
    os.path.expanduser('~'), 'femag')
try:
    os.makedirs(workdir)
except OSError:
    pass

femag = femagtools.Femag(workdir)

operatingConditions = dict(
    angl_i_up=0.0,
    calculationMode="torq_calc",
    wind_temp=20.0,
    magn_temp=20.0,
    current=7.07,
    poc=poc,
    speed=10.0)

r = femag(machine,
          operatingConditions)

print("""
Force [N] = {}
""".format(r.linearForce_fft[0]['force'][0]))
print("""
Losses [W]:

 Stator Teeth {0:8.1f}
        Yoke  {1:8.1f}
 Rotor        {2:8.1f}
 Magnet       {3:8.1f}
 Windings     {4:8.1f}

        Total {5:8.1f}
""".format(
    r.losses[-1]['staza'],
    r.losses[-1]['stajo'],
    r.losses[-1]['rotfe'],
    r.losses[-1]['magnetJ'],
    r.losses[-1]['winding'],
    r.losses[-1]['total']))
