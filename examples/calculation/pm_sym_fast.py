import femagtools
import femagtools.poc
import pathlib
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
        num_layers=1,
        leak_dist_wind=dict(
            perimrad=67.1e-3,  # Radius of perimeter [m]
            vbendrad=5e-3,    # Bending radius vertical [m]
            endheight=20e-3,  # End winding height [m]
            wiredia=1e-3)     # Wire diameter [m]
    )
)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

workdir = pathlib.Path().home() / 'femag'
workdir.mkdir(parents=True, exist_ok=True)

femag = femagtools.Femag(workdir)

poc = femagtools.poc.HspPoc(harm=[1, 5],
                            amp=[1, 0.01],
                            phi=[0, 0])

simulation = dict(
    angl_i_up=0.0,
    calculationMode="pm_sym_fast",
    wind_temp=60.0,
    magn_temp=60.0,
    current=50.0,
    speed=50.0,
    period_frac=6,
    poc=poc,
    plots=['field_lines', ['Babs', 0.5, 3.2, 'Babs.svg']])

r = femag(machine,
          simulation)

print('Torque [Nm] = {}'.format(r.machine['torque']))
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
