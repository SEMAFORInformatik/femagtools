import femagtools
import femagtools.poc
import femagtools.amela
import pathlib
import logging
import math

machine = dict(
    name="PM130L4",
    lfe=30e-3,
    poles=10,
    outer_diam=0.1,
    bore_diam=0.055,
    inner_diam=0.035,
    airgap=0.001,
    stator=dict(
        num_slots=12,
        statorfsl=dict(
            sw=0.005,
            tw=0.007,
            slot_h1=0.0017,
            slot_h2=0.002)
    ),
    magnet=dict(
        magnetfsl=dict(
           hm=3e-3,
           bm=11.205e-3
        )
    ),

    winding=dict(
        num_phases=3,
        num_wires=10,
        coil_span=1,
        num_layers=2)
)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

workdir = pathlib.Path().home() / 'femag'
workdir.mkdir(parents=True, exist_ok=True)

femag = femagtools.Femag(workdir)

simulation = dict(
    angl_i_up=0.0,
    calculationMode="pm_sym_fast",
    wind_temp=60.0,
    magn_temp=20.0,
    current=30.0/math.sqrt(2),
    speed=8000/60,
    num_move_steps=-73,
    period_frac=1,
    plots=['field_lines', ['Babs', 0.5, 3.2, 'Babs.svg']])

r = femag(machine,
          simulation)

# amela directory can be specified separately
amela_dir = '../../amela_interface'
amela = femagtools.amela.Amela(workdir, dict(name=machine['name']), amela_dir)
# magnet_loss is a python dict that contains the loss data
magnet_loss = amela()
