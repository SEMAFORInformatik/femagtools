# -*- coding: utf-8 -*-
import femagtools
import matplotlib.pylab as pl
import femagtools.plot

bch = femagtools.read_bchfile('TEST_001.BCH')

beta = bch.ldq['beta']
i1 = bch.ldq['i1']
torque = bch.ldq['torque']

femagtools.plot.i1beta_torque_plot(i1, beta, torque)
pl.show()

