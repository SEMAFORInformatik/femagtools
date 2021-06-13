# Run AC Simulation
# Ronald Tanner Semafor Informatik & Energie AG Basel
import logging
import json
import femagtools

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

machine = dict(
    name="VTT-PP",
    lfe=0.16,
    poles=4,
    outer_diam=0.220,
    bore_diam=0.125,
    inner_diam=0.056,
    airgap=0.5e-3,

    stator=dict(
        num_slots=48,
        mcvkey_yoke="M800-65A",
        nodedist=2.0,
        statorRotor3=dict(
            slot_height=0.0283,
            slot_h1=7e-4,
            slot_h2=9e-4,
            slot_r1=0.0,
            slot_r2=3.4e-3,
            wedge_width1=4e-3,
            wedge_width2=0,
            tooth_width=0.0,
            slot_top_sh=0.0,
            slot_width=2.8e-3
        )
    ),

    rotor=dict(
        num_slots=40,
        mcvkey_yoke="M800-65A",
        sigma=24e6,  # Al conductivity S/m
        statorRotor3=dict(
            slot_height=0.0156,
            slot_h1=0.5e-3,
            slot_h2=2.2e-3,
            slot_r1=2.2e-3,
            slot_r2=1e-3,
            wedge_width1=1.0e-3,
            wedge_width2=0,
            tooth_width=0.0,
            slot_top_sh=0.0,
            slot_width=1e-3,
            conductivity=24e6  # Al
        )

    ),

    windings=dict(
        dia_wire=1.5e-3,
        num_phases=3,
        num_layers=1,
        num_wires=32,
        coil_span=12,
        fillfac=0.42,  #
        rellength=1.4  #
    )
)

workdir = '.'
femag = femagtools.Femag(workdir, magnetizingCurves='.')
#                         cmd='xfemag64.sh')

# set TS simulation
simulation = dict(
    calculationMode="asyn_motor",
    wind_temp=60,
    bar_temp=120,
    speed=1450/60,
    f1=50.0,
    num_par_wdgs=2,
    wdgcon=1,  # 0:open, 1:star, 2:delta
    u1=230.0)  # phase voltage

r = femag(machine, simulation)

print(json.dumps(r))
