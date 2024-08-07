#!/usr/bin/env python
# -*- coding: utf-8 -*-
import femagtools
import femagtools.machine
import numpy as np
import matplotlib.pylab as pl
import warnings
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

warnings.simplefilter("error")
bch = femagtools.read_bchfile('../bch-erg/TEST_001.BCH')

beta = bch.ldq['beta']
i1 = bch.ldq['i1']
torque = bch.ldq['torque']

r10 = 0.1
p = 4

pm = femagtools.machine.create(bch, r1=r10, ls=1e-3)

u1 = 340
tq = 200

iqx, idx = pm.iqd_torque(tq)
w1 = pm.w1_u(u1, iqx, idx)
i1 = femagtools.machine.betai1(iqx, idx)[1]

fig, ax = pl.subplots()

id = np.linspace(pm.iqdmin(i1)[1],
                 pm.iqdmax(i1)[1])

iq = [pm.iq_u(w1, u1, ix) for ix in id]
ax.plot(id, iq, label='U1={} V'.format(u1))

i1x = pm.i1_torque(tq, beta[-1] / 180 * np.pi)
iq = np.linspace(0.8*iqx, 1.8*iqx, 20)
id = np.array([pm.id_torque(tq, ix) for ix in iq])
ax.plot(id, iq, label='Tq={} Nm'.format(tq))

ax.annotate('f1={0:4.1f} Hz'.format(w1 / np.pi / 2),
            xy=(idx, iqx), xytext=(1.3 * idx, 1.5 * iqx),
            arrowprops=dict(arrowstyle="->"), fontsize=12)

ax.arrow(0, 0, idx + 0.075 * i1, iqx - 0.08 * i1, color='r',
         head_width=0.05 * i1, head_length=0.08 * i1)
ax.text(1.38 * idx, 0.5 * iqx,
        r'$I_1={0:3.1f} A$'.format(
            femagtools.machine.betai1(iqx, idx)[1]), fontsize=14)
ax.arrow(0, 0, 0, 170, color='k', head_width=0.05 * i1, head_length=0.08 * i1)

ax.annotate("",
            xy=(0.36 * idx, 0.36 * iqx), xycoords='data',
            xytext=(0, 0.4*i1), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            color="0.5",
                            shrinkA=5, shrinkB=5,
                            patchA=None,
                            patchB=None,
                            connectionstyle="arc3,rad=0.3"))
ax.text(0.52 * idx, 0.58 * iqx,
        r'$\beta={0:3.1f}^o$'.format(np.arctan2(idx, iqx)/np.pi * 180),
        fontsize=12)

ax.grid()
xlim = ax.get_xlim()
ax.set_xlim([xlim[0], 20])
ax.set_xlabel('Id / A')
ylim = ax.get_ylim()
ax.set_ylim([0, ylim[1]])
ax.set_ylabel('Iq / A')
legend = ax.legend(loc='upper center')

ax.set_aspect('equal')

pl.show()
