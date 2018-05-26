#!/usr/bin/env python
# -*- coding: utf-8 -*-
import femagtools.forcedens
import matplotlib.pylab as pl

fdens = femagtools.forcedens.read('example.PLT1')

f, (ax1, ax2) = pl.subplots(2, sharex=True)
ax1.set_title('{}, Rotor position {}'.format(
    fdens.title, fdens.positions[0]['position']))
ax1.plot(fdens.positions[0]['X'], [1e-3*ft
                                   for ft in fdens.positions[0]['FT']],
         label='F tang')
ax1.plot(fdens.positions[0]['X'], [1e-3*ft
                                   for ft in fdens.positions[0]['FN']],
         label='F norm')
ax1.legend()
ax2.plot(fdens.positions[0]['X'], fdens.positions[0]['B_T'],
         label='B tang')
ax2.plot(fdens.positions[0]['X'], fdens.positions[0]['B_N'],
         label='B norm')
ax1.set_ylabel('kN/m²')
ax2.set_ylabel('T')
pl.xlabel('Position / deg')
ax2.legend()
pl.show()
