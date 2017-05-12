# -*- coding: utf-8 -*-
import femagtools
import logging

mcvData = dict(curve=[dict(
    bi=[0.0, 0.09, 0.179, 0.267, 0.358,
        0.45, 0.543, 0.6334, 0.727,
        0.819, 0.9142, 1.0142, 1.102,
        1.196, 1.314, 1.3845, 1.433,
        1.576, 1.677, 1.745, 1.787,
        1.81, 1.825, 1.836],
        
    hi=[0.0, 22.16, 31.07, 37.25, 43.174,
        49.54, 56.96, 66.11, 78.291,
        95, 120.64, 164.6, 259.36,
        565.86, 1650.26, 3631.12, 5000, 10000,
        15000, 20000, 25000, 30000, 35000, 40000])],
               desc=u"Demo Steel",
               ch=4.0,
               cw_freq=2.0,
               cw=1.68)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

mcv = femagtools.mcv.MagnetizingCurve(mcvData)
mcv.writefile('m270-35a')
