# FEMAG DC PM/Rel-Simulation
#
import logging
import femagtools
import pathlib
from machine import machine
import matplotlib.pyplot as plt

simulation = dict(
    calculationMode="calc_field_ts",
    src_ampl=10.0,
    src_type='current',
    freq=100.0,
    speed=50.0,
    sim_time=[0.01, 0.02],
    store_time=0.01/20,
    dtmin=0.01/20,
    dtmax=0.01/20,
    resmin=1e-5,
    resmax=1e-4,
    vtu_movie=True)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')
modelname = machine['name'].replace(' ', '_')
workdir = pathlib.Path.home() / 'femag'

femag = femagtools.Femag(workdir, magnetizingCurves='../magnetcurves',
                         cmd='tsfemag64',
                         templatedirs=['templates'])

r = femag(machine, simulation)

fig, axs = plt.subplots(2, 2, sharex=True)
axs[0, 0].plot(r['time'], r['I1'])
axs[0, 0].plot(r['time'], r['I2'])
axs[0, 0].plot(r['time'], r['I3'])
axs[0, 0].set_title('Current / A')
axs[0, 0].grid()
axs[0, 1].plot(r['time'], r['U1'])
axs[0, 1].plot(r['time'], r['U2'])
axs[0, 1].plot(r['time'], r['U3'])
axs[0, 1].set_title('Voltage / V')
axs[0, 1].grid()
axs[1, 0].plot(r['time'], r['angle'])
axs[1, 0].set_title('Angle / °')
axs[1, 0].grid()
axs[1, 0].set_xlabel('Time / s')
axs[1, 1].plot(r['time'], r['torque'])
axs[1, 1].set_title('Torque / Nm')
axs[1, 1].grid()
axs[1, 1].set_xlabel('Time / s')
fig.tight_layout()
plt.show()
