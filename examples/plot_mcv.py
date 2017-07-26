#!/usr/bin/env python
# read mc file with filename given as command line argument and
# sho BH and Mue plot
#
# Usage:
#  plot.mcv.py <filename>
#
import matplotlib.pylab as pl
import femagtools.mcv
import math
import sys

MUE0 = 4e-7*math.pi

mcv = femagtools.mcv.Reader()
mcv.readMcv(sys.argv[1])
r = mcv.get_results()

bh = [(bi, hi)
      for bi, hi in zip(r['curve'][0]['bi'],
                        r['curve'][0]['hi']) if bi > 0 and hi > 0]

ji = [b-MUE0*h for b, h in bh]
bi, hi = zip(*bh)

fig = pl.figure()
g = pl.GridSpec(1, 2)
ax = [pl.subplot(g[0, 0]),
      pl.subplot(g[0, 1])]
ax[0].set_title(r['name'])
ax[0].semilogx(hi, bi, label='Induction')
ax[0].semilogx(hi, ji, label='Polarisation')
ax[0].set_xlabel('H / A/m')
ax[0].set_ylabel('T')
ax[0].legend(loc='lower right')
ax[0].grid()
    
ur = [bx/hx/MUE0 for bx, hx in bh]
ax[1].plot(bi, ur)
ax[1].set_xlabel('B / T')
ax[1].set_title('Muer')
ax[1].grid()
    
fig.tight_layout()
fig.subplots_adjust(top=0.94)
pl.show()

