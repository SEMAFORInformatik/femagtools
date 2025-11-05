# FEMAG TS 3ph short circuit
#
import logging
import femagtools
from femagtools import netlist
import pathlib
from machine import machine
import matplotlib.pyplot as plt

speed = 50
freq = speed*machine['poles']/2
cond_temp=140
current_angles=[60,180,300]

current = 100  # amplitude of phase current in A
simulation = dict(
    calculationMode="calc_field_ts",
    speed=speed,
    magn_temp=80,
    cond_temp=120,
    netlists=[
        netlist.create_shortcircuit(
            freq, current, current_angles, [0, 0],
            phseq='123'),
        netlist.create_shortcircuit(
            freq, 0, current_angles, [1, 1], phseq='123')],
    sim_time=[0.01,0.02],
    dtmin=0.01/40,
    dtmax=0.01/40,
    resmin=1e-5,
    resmax=1e-4,
    vtu_movie=True)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')
modelname = machine['name'].replace(' ', '_')
workdir = pathlib.Path.home() / 'femag'

femag = femagtools.Femag(workdir, magnetizingCurves='../magnetcurves',
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
