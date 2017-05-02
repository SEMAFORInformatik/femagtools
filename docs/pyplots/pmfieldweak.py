#!/usr/bin/env python
# -*- coding: utf-8 -*-
import femagtools
import femagtools.machine
import numpy as np
import matplotlib.pylab as pl

ld = 0.0014522728
lq = 0.0038278836
psim = 0.11171972
p = 4

pm = femagtools.machine.PmRelMachineLdq(3, p,
                                        ld=ld,
                                        lq=lq,
                                        psim=psim)

u1 = 340
tq = 170

iqx, idx = pm.iqd_torque(tq)
w1 = pm.w1_u(u1, iqx, idx)

fig, ax = pl.subplots()

id = np.linspace(-120, -20)
iq = [pm.iq_u(w1, u1, ix) for ix in id]
ax.plot(id, iq, label='U1={} V'.format(u1))

iq = np.linspace(65, 180)
id = np.array([pm.id_torque(tq, ix) for ix in iq])
ax.plot(id, iq, label='Tq={} Nm'.format(tq))

i1 = np.linalg.norm(np.array((iqx, idx)))/np.sqrt(2)
ax.annotate('f1={0:4.1f} Hz'.format(w1 / np.pi / 2),
            xy=(idx, iqx), xytext=(1.3 * idx, 1.5 * iqx),
            arrowprops=dict(arrowstyle="->"))

ax.arrow(0, 0, idx + 0.075 * i1, iqx - 0.08 * i1, color='r',
         head_width=0.05 * i1, head_length=0.08 * i1)
ax.text(1.38 * idx, 0.5 * iqx,
        r'$I_1={0:3.1f} A$'.format(i1), fontsize=18)
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
        fontsize=14)

ax.grid()
xlim = ax.get_xlim()
ax.set_xlim([xlim[0], 10])
ax.set_xlabel('Id / A')
ylim = ax.get_ylim()
ax.set_ylim([0, ylim[1]])
ax.set_ylabel('Iq / A')
legend = ax.legend(loc='upper center')

ax.set_aspect('equal')

pl.show()
