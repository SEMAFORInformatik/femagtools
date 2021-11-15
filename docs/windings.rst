Windings
********

The windings module provides utility functions and the class
Windings to define and analyze windings.

Example of a 3-phase, 2-layer Winding with 54 slots and 6 pole pairs::

  w = femagtools.windings.Winding(dict(Q=54, p=6, m=3, l=2))

  assert wdg.yd == 4
    assert wdg.zoneplan() == (
        [[1, 2, -6], [4, 5, -9], [-3, 7, 8]],
        [[1, -5, -6], [4, -8, -9], [-2, -3, 7]])

  femagtools.plot.zoneplan(w)
  
.. figure:: img/zoneplan.png

  femagtools.plot.winding(w)
  
.. figure:: img/winding.png

  mmf = w.mmf()
  femagtools.plot.mmf(mmf)
  femagtools.plot.mmf_fft(mmf)

.. figure:: img/mmf.png
            
.. figure:: img/mmf_fft.png
