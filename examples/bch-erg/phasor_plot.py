import femagtools
import femagtools.plot
import matplotlib.pylab as pl

r = femagtools.read_bchfile('TEST_002.BCH')

femagtools.plot.phasor(r)
pl.show()

#w1 = r.dqPar['npoles']*r.dqPar['speed']*math.pi
#femagtools.plot.phasor(r.dqPar['up'], r.dqPar['i1'][-1],
#                       r.dqPar['beta'][-1], 0,
#                       w1*r.dqPar['ld'][-1],
#                       w1*r.dqPar['lq'][-1])
