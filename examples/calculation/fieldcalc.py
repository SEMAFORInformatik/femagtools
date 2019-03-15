import femagtools
import femagtools.airgap
import os
import matplotlib.pyplot as plt
import logging


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


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
    
    windings=dict(
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

simulation = dict(
    calculationMode="fieldcalc")

r = femag(machine,
          simulation)

# read airgap flux density 
Q = machine['stator']['num_slots']
p = machine['poles']//2
taup = 360/gcd(Q, p)

rag = femagtools.airgap.read(
    os.path.join(workdir, 'bag.dat'), taup)

plt.plot(rag['pos'], rag['B'])
plt.grid()
plt.show()
