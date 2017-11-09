#!/usr/bin/env python
# read mc file with filename given as command line argument and
# sho BH and Mue plot
#
# Usage:
#  plot.mcv.py <filename>
#
import matplotlib.pylab as pl
import femagtools.mcv
import femagtools.plot
import sys

mcv = femagtools.mcv.Reader()
mcv.readMcv(sys.argv[1])

if mcv['mc1_type'] in (femagtools.mcv.MAGCRV, femagtools.mcv.ORIENT_CRV):
    ncols = 2
else:  # Permanent Magnet
    ncols = 1

fig, ax = pl.subplots(nrows=1, ncols=ncols)
if ncols > 1:
    pl.subplot(1, ncols, 1)
    femagtools.plot.mcv_hbj(mcv)
    pl.subplot(1, ncols, 2)
    femagtools.plot.mcv_muer(mcv)
else:
    pl.subplot(1, ncols, 1)
    femagtools.plot.mcv_hbj(mcv, log=False)
    
fig.tight_layout()
fig.subplots_adjust(top=0.94)
pl.show()
