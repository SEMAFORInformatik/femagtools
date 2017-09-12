#!/usr/bin/env python
# -*- coding: utf-8 -*-
import femagtools
import matplotlib.pylab as pl
import numpy as np
import scipy.interpolate as ip

bch = femagtools.read_bchfile('TEST_002.BCH')

pos = bch.torque[-1]['angle']
torque = bch.torque[-1]['torque']

fig = pl.figure()
ax = fig.add_subplot(111)
ax.set_title(bch.project)
ax.set_xlabel(u'Position [°]')
ax.set_ylabel(u'Torque [Nm]')
ax.plot(pos, torque, 'go', label='Calculated')

k = 20
pos = np.linspace(0, max(bch.torque[-1]['angle']),
                  k*len(bch.torque[-1]['torque']))
f = ip.interp1d(bch.torque[-1]['angle'],
                bch.torque[-1]['torque'], kind='cubic')

ax.plot(pos, f(pos), label='Interpolated')
ax.plot(pos, [np.average(
    bch.torque[-1]['torque'])]*len(pos), label='Average')
ax.legend()

pl.grid()
pl.show()
