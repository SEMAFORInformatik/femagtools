import femagtools
import femagtools.plot.bch
import femagtools.machine
import femagtools.shortcircuit
import logging
import matplotlib.pyplot as plt
import os

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

machine = dict(
    name="PM 270 L8",
    lfe=0.08356,
    poles=8,
    outer_diam=0.26924,
    bore_diam=0.16192,
    inner_diam=0.092,
    airgap=0.00075,

    stator=dict(
        num_slots=48,
        nodedist=2.5,
        mcvkey_yoke='M330-50A',
        statorRotor3=dict(
            slot_height=0.0335,
            slot_h1=0.001,
            slot_h2=0.0,
            slot_r1=0.0001,
            slot_r2=0.00282,
            wedge_width1=0.00295,
            wedge_width2=0,
            middle_line=0,
            tooth_width=0.0,
            slot_top_sh=0.0,
            slot_width=0.00193)
    ),

    magnet=dict(
        mcvkey_yoke='M330-50A',
        magnetIronV=dict(
            magn_width=18e-3,
            magn_height=6.48e-3,
            magn_angle=145,
            magn_num=1,
            gap_ma_iron=0.2e-3,
            air_triangle=1e-3,
            iron_height=2.61e-3,
            iron_hs=0.1e-3,
            shaft_rad=55.32e-3,
            iron_shape=80.2e-3,
            air_space_h=5.5e-3,
            iron_bfe=3e-3,
            magn_di_ra=6e-3,
            corner_r=0,
            air_sp_ori=1,
            magn_ori=1,
            condshaft_r=55.32e-3)
    ),

    winding=dict(
        num_phases=3,
        num_wires=9,
        coil_span=6.0,
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

femag = femagtools.Femag(workdir, magnetizingCurves='../magnetcurves')

pmRelSim = dict(
    angl_i_up=-39.3,
    calculationMode="pm_sym_fast",
    wind_temp=60.0,
    magn_temp=60.0,
    current=76.43,
    period_frac=6,
    speed=50.0)

r = femag(machine, pmRelSim)
scsim = dict(
    l_end_winding=0,
    l_external=0,
    sc_type=3,
    initial=2,
    allow_demagn=0,
    sim_demagn=1)
sc3 = femagtools.shortcircuit.shortcircuit(
    femag, machine, r, scsim)
print('Torque [Nm] = {}'.format(r.machine['torque']))
print('''
Short Circuit    Current         Torque
 Peak       iks {2:8.1f} A  tks {3:8.1f} Nm
 Stationary ikd {0:8.1f} A  tkd {1:8.1f} Nm

  peak winding currents {4}
'''.format(sc3['ikd'],
           sc3['tkd'],
           sc3['iks'],
           sc3['tks'],
           sc3['peakWindingCurrents']))
print(f'Demag {sc3["demag"]}')
fig, axs = plt.subplots(nrows=4, figsize=(10, 12), sharex=True)
femagtools.plot.bch.transientsc_currents(sc3, ax=axs[0],
                                         title='3 phase fault',
                                         set_xlabel=False)
femagtools.plot.bch.transientsc_torque(sc3, ax=axs[1],
                                       set_xlabel=False)
from femagtools.multiproc import Engine
engine = Engine()
scsim['sc_type'] = 2
sc2 = femagtools.shortcircuit.shortcircuit(
    femag, machine, r, scsim, engine)

print('''
Short Circuit    Current         Torque
 Peak       iks {2:8.1f} A  tks {3:8.1f} Nm
 Stationary ikd {0:8.1f} A  tkd {1:8.1f} Nm

  peak winding currents {4}
'''.format(sc2['ikd'],
           sc2['tkd'],
           sc2['iks'],
           sc2['tks'],
           sc2['peakWindingCurrents']))
print('Demag {}'.format(sc2['demag']))

femagtools.plot.bch.transientsc_currents(sc2, ax=axs[2],
                                         title='2 phase fault',
                                         set_xlabel=False)
femagtools.plot.bch.transientsc_torque(sc2, ax=axs[3])
fig.tight_layout()
fig.savefig('shortcircuits.pdf')

fig, axs = plt.subplots(nrows=2)
femagtools.plot.bch.demagnetization(sc3['demag'], ax=axs[0])
femagtools.plot.bch.demagnetization(sc2['demag'], ax=axs[1])
fig.tight_layout()
fig.savefig('demagnetization.pdf')
