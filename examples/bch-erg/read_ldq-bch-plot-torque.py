# -*- coding: utf-8 -*-
import femagtools
import matplotlib.pylab as pl
import femagtools.plot
import femagtools.machine

bch = femagtools.read_bchfile('TEST_001.BCH')

beta = bch.ldq['beta']
i1 = bch.ldq['i1']
torque = bch.ldq['torque']

pm = femagtools.machine.create(bch, r1=0, ls=0)
femagtools.plot.mtpa(pm, i1[-1])
pl.show()

