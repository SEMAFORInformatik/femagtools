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
posip = np.linspace(0, pos[-1], k*len(torque))

f = ip.interp1d(pos, torque, kind='cubic')

ax.plot(posip, f(posip), label='Interpolated')
ax.plot(pos, [np.average(torque)] * len(pos), label='Average')
ax.legend()

pl.grid()
pl.show()
