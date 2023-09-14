.. _sizing:

**Sizing**
**********

The machine sizing module includes functions that
creates a machine model dict (see :ref:`machine-model`) for following types of
radial flux machines:

* Surface Mounted PM (spm)
* Interior Magnet PM (ipm)
* Electric Excited SM (eesm)
* Induction Machine (im)

Each function has three mandatory parameters:

* pnom: nominal (rated) power in W
* speed: shaft speed in 1/s
* p: pole pairs

In addition a couple of optional parameters with
default values are available:

==============  ======================================  ==========  ======
Name            Description                             Unit        Value
==============  ======================================  ==========  ======
airgap          airgap width                            m           1.5e-3
eta             efficiency                                          0.92
cos_phi         power factor                                        0.95
m               number of phases                                    3
ui_u            ratio U ind / U                                     0.62
lda             length/taup ratio                                   1.33
J               current density 3 .. 6                  A/m²        3.8e6
sigmas          shear force 10 .. 45 kN/m2              N/m²        17e3
Ba              flux density in airgap                  T           0.7
Bth             flux density in teeth                   T           1.5
By              flux density in yoke                    T           1.2
Bthr            flux density in rotor teeth             T           1.7
Byr             flux density in rotor yoke              T           1.4
kq              stator winding fill factor                          0.42
kqr             rotor winding fill factor                           0.6
mag_width       rel magnet width                                    0.8
Hc              max. coercitive field strength          kA/m        700
brem            remanence                               T           1.15
demag           safety factor for demagnetisation                   1.5
external_rotor  external rotor                                      False
coil_span       coil span                                           0
hfe             iron height between magnet and airgap   m           1e-3
==============  ======================================  ==========  ======

For IM the default values are as following:

==============  ======================================  ==========  ========
Name            Description                             Unit        Value
==============  ======================================  ==========  ========
airgap          airgap width                            m           0.75e-3
eta             efficiency                                          0.87
cos_phi         power factor                                        0.8
==============  ======================================  ==========  ========


Sizing Example for a SPM::

  p2 = 1.5e3
  speed = 1500/60
  udc = 550
  p = 4

  machine = femagtools.machine.sizing.spm(
     p2, speed, p, udc=udc, Q1=36,
     Hc=700, sigmas=12e3, brem=1.1, Ba=0.77,
     cos_phi=0.7, eta=0.8, demag=1.7, lda=0.9)
