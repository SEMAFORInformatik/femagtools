
import femagtools.erg
import femagtools.machine
import femagtools.plot

import matplotlib.pylab as pl

res = femagtools.erg.read('ldlq.erg')

lfe = 350

pm = femagtools.machine.create(res, r1=0, ls=0, lfe=lfe)

femagtools.plot.mtpa(pm, res['i1'][-1])
pl.show()
