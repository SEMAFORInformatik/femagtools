MagnetizingCurve
****************

The MagnetizingCurve is a container of magnetizing curves (eg. lamination material) that can be referenced by the model mcvkey attributes. It can ether point to a directory of MC/MCV-File or hold a list of magnet curves which are identified by name.

Each magnetizing curve is described by the following properties

================  ================================ ======== =======
Attribute          Description                     Unit     Default
================  ================================ ======== =======
   name           Identifier of this curve
   desc           Description
   curve          List of dictionaries with
                  bi (list of induction values)    T,
                  hi (List of field strength       A/m
		  values)                          
   ch             hysteresis loss factor                    0
   cw             eddy current loss factor                  0
   ch_freq        hysteresis exponent                       0
   cw_freq        eddy-current exponent                     0
   b_coeff        induction loss exponent                   0
   Bo             reference induction              T        1.5
   fo             reference frequency              Hz       50
   fillfac        iron fill factor                          1
   bsat           saturation induction             T        2.15
   rho            specific weight                  kg/dm3   7.65
================  ================================ ======== =======

Loss calculation formula:

 (cw*(f/fo)**cw_freq + ch*(f/fo)**ch_freq)*(B/Bo)**b_coeff

The Reader object which is included in the mcv module can be used to read MCV/MC files.

Permeability and polarisation calculation example::

  MUE0 = 4e-7*math.pi

  mcv = femagtools.mcv.Reader()
  mcv.readMcv('magnetcurves/M270-35A.MCV')
  r = mcv.get_results()

  bh = [(bi, hi)
        for bi, hi in zip(r['curve'][0]['bi'],
                          r['curve'][0]['hi']) if bi > 0 and hi > 0]

  ji = [b-MUE0*h for b, h in bh]
  muer = [bx/hx/MUE0 for bx, hx in bh]


Using a magnetizingcurve to Write a mcv file::

   mcvData = dict(curve=[ dict(
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
           15000, 20000, 25000, 30000, 35000, 40000]
       )],
       desc=u"Demo Steel",
       ch=4.0
       cw_freq=2.0,
       cw=1.68)

    mcv = femagtools.mcv.MagnetizingCurve(mcvData)
    
    mcv.writefile('m270-35a')

.. image:: img/mcv.png
  :height: 290pt
