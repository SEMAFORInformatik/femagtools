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

if r['ctype'] in (femagtools.mcv.MAGCRV, femagtools.mcv.ORIENT_CRV):
    ncols = 2
    bh = [(bi, hi*1e-3)
          for bi, hi in zip(r['curve'][0]['bi'],
                            r['curve'][0]['hi']) if bi > 0 and hi > 0]
    
    ji = [b-MUE0*h*1e3 for b, h in bh]
else:  # Permanent Magnet
    ncols = 1
    bh = [(bi, hi*1e-3)
          for bi, hi in zip(r['curve'][0]['bi'],
                            r['curve'][0]['hi'])]
    ji = []
bi, hi = zip(*bh)

fig, ax = pl.subplots(nrows=1, ncols=ncols)
if ncols > 1:
    ax[0].set_title(r['name'])
    ax[0].semilogx(hi, bi, label='Induction')
    ax[0].semilogx(hi, ji, label='Polarisation')
    ax[0].set_xlabel('H / kA/m')
    ax[0].set_ylabel('T')
    ax[0].legend(loc='lower right')
    ax[0].grid()

    ur = [bx/hx/MUE0*1e-3 for bx, hx in bh]
    ax[1].plot(bi, ur)
    ax[1].set_xlabel('B / T')
    ax[1].set_title('rel. Permeability')
    ax[1].grid()
else:
    ax.set_title(r['name'])
    ax.plot(hi, bi, label='Induction')
    ax.set_xlabel('H / kA/m')
    ax.set_ylabel('T')
    ax.grid()
    
fig.tight_layout()
fig.subplots_adjust(top=0.94)
pl.show()
