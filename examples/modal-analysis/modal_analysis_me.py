import femagtools
import femagtools.plot
import matplotlib.pyplot as plt
import pathlib
import logging

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
        num_slots_gen=12,
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

    windings=dict(
        num_phases=3,
        num_wires=10,
        coil_span=1,
        num_layers=2)
)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

workdir = pathlib.Path().home() / 'femag'
workdir.mkdir(parents=True, exist_ok=True)

# remove old results (if any)
for i in workdir.glob('*_Mode_*'): 
    i.unlink()

femag = femagtools.Femag(workdir)

# mechanical material can be passed in a list: 
# [mass_density, Young's Modulus, Possion's number]
#  kg/m3, Gpa, #
simulation = dict(
    num_modes=15,
    calculationMode="modal_analysis",
    stator_material=[7700, 210, 0.3], 
    slot_material=[5000, 1.5, 0.3],
    export_figure=True
)
# return eigenfrequency, eigenvectors
r = femag(machine,
          simulation)
          
if simulation['export_figure']:
    femagtools.plot.eigenmode(r)
    plt.show()



