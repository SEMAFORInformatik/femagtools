BchReader
*********

The BchReader object holds the most important FEMAG results. It has
following attributes which mostly correspond to the text sections in the file:

================  =======================================================
Attribute          Description     
================  =======================================================
project            Name of model file
filename           Name of BCH file
date               calculation date
version            FEMAG version
nodes              number of nodes
elements           number of elements
quality            meshing quality
flux               Flux observed
flux_fft           Fourier-Analysis of flux values
torque             Torque-Force values
torque_fft         Fourier-Analysis of torque values
linearForce        Force x and y values
linearForce_fft    Fourier-Analysis of force values
psidq              PSID-Psiq-Identification
psidq_ldq          PSID-Psiq-Identification (Ld, Lq)
machine            Machine data
lossPar            Control parameters for Loss calculation
magnet             Magnet data
airgapInduction    airgap induction
scData             Transient short circuit
dqPar              DQ-Parameter for open Winding Modell
ldq                Ld-Lq-Identification
losses             Losses in iron, magnets and conductors from Move-calc
demag              Demagnetisation
weights            Total weight and weight of iron, conductor and magnets
================  =======================================================

Flux
====

* Flux: list of dictionaries for each winding

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  displ             position
  displunit         unit (mm, deg) of position values
  flux_k            flux 
  voltage_dpsi      voltage dpsi/dt
  voltage_four      voltage (fourier transformation)
  current_k         current
  voltage_ir        voltage
  ================  =======================================================

  
* flux_fft: list of dictionaries for each winding

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  order             order of harmonic
  flux              flux amplitude
  flux_perc         flux amplitude percentage of base harmonic
  voltage           voltage amplitude
  voltage_perc      voltage amplitude percentage of base harmonic
  a                 amplitude of sin term
  b                 amplitide of cos term
  ================  =======================================================

Torque
======

* torque: list of dictionaries for each current and/or beta angle

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  angle             Position
  current_1         Current
  force_x           Force in tangetial direction 
  force_y           Force in radial direction 
  t_idpsi           Torque
  torque            Torque
  ripple            Diff between max and min torque value
  ================  =======================================================


* torque_fft: list of dictionaries for each current and/or beta angle

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  order             order of harmonic list
  torque            Torque amplitude list
  torque_perc       Torque in percentage of first amplitude
  a                 amplitude list of sin term
  b                 amplitude list of cos term
  ================  =======================================================

Linear Force
============

* linearForce: list of dictionaries for each current and/or beta angle

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  displ             Position
  magnet_1          Current
  force_x           Force in x direction 
  force_y           Force in y direction 
  ================  =======================================================


* linearForce_fft: list of dictionaries for each current and/or beta angle in x- and y-direction

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  order             order of harmonic list
  force             Force amplitude list
  force_perc        Force in percentage of first amplitude
  a                 amplitude list of sin term
  b                 amplitude list of cos term
  ================  =======================================================

Psidq
=====

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  iq                Iq current list (n)
  id                Id current list (m)
  psid              Psid matrix (n x m)
  psiq              Psiq matrix (n x m)
  torque            Torque matrix (n x m)
  ================  =======================================================
  
Psidq Ldq
=========

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  iq                Iq current list (n)
  id                Id current list (m)
  ld                Ld matrix (n x m)
  lq                Lq matrix (n x m)
  psim              Psim matrix (n x m)
  psid              Psid matrix (n x m)
  psiq              Psiq matrix (n x m)
  torque            Torque matrix (n x m)
  ================  =======================================================

Machine
=======

  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  beta              Beta list                                   deg
  plfe1             Iron losses stator                          W
  plfe2             Iron Losses rotor                           W
  plmag             Magnet losses                               W
  plcu              Winding losses                              W
  pltotal           Total losses                                W
  plfe              Total Iron losses                           W
  lfe               Length of armature                          m
  eff               Efficiency                                  %
  m                 Number of phases
  p                 Number of pole pairs
  Q                 Number of stator slots
  p2                Mechanical power                            W
  i1                Phase current                               A
  A                 current loading                             kA/m
  J                 current density                             A/mm2
  kcu               copper fill factor                          %
  AJ                Therm loading                               A/cm.mm2
  torque            Torque                                      Nm
  fd                Force density                               N/mm²
  ld                Ld Inductance                               H
  lq                Lq Inductance                               H
  r1                Stator resistance                           Ohm
  psim              Magn flux                                   Vs
  n                 Speed                                       1/s
  lpfe1_0           Iron Losses in stator at noload             W
  lpfe2_0           Iron Losses in rotor at noload              W
  lpmag_0           Magnet losses at noload                     W
  ================  ========================================== =============
  
  Example::
    
    {'m': 3,
    'p': 4,
    'qs_sim': 12,
    'p_sim': 2,
    'Q': 48,
    'n': 50.0,
    
    'kcu': 40.0,
    'r1': 0.055,
    'AJ': 84365.4609,
    'A': 213.2994,
    'fd': 119.0008,
    'J': 39.5526,
    
    'lfe': 0.08356,
    'ld': 0.0008625,
    'lq': 0.00132,
    'psim': 0.1152,

    'torque': 405.7295,
    'p2': 127463.7,

    'plfe1_0': 172.9209,
    'plmag_0': 0.0239,
    'plfe2_0': 0.7076,
    'i1': 500.0,
    'beta': [0.0, -25.0],

    'plfe1': [1463.3809, 1374.8728],
    'plfe2': [71.727, 77.0296],
    'plmag': [4.1524, 15.1965],
    'plcu': [10305.4824, 10305.4824],
    'pltotal': [11844.7427, 11772.581300000002],
    'plfe': [1535.1079000000002, 1451.9024000000002]
    'eff': 91.5449}

DqPar
=====

  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  beta              Beta list                                   deg
  lfe               Length of armature                          m
  npoles            Number of poles
  cosphi            Power factor
  ld                Inductance Ld                               H
  lq                Inductance Lq                               H
  psid              Flux in d-axis                              Vs
  psiq              Flux in q-axis                              Vs
  psim              Magnetizing Flux                            Vs
  psim0             Magnetizing Flux at no-load                 Vs
  u1                Terminal voltage                            V
  up                MMF voltage                                 V
  up0               MMF voltage at-noload                       V
  u1                Terminal voltage                            V
  gamma             Angle between Up and U1                     deg
  i1                Phase current                               A
  phi               Angle between U1 and I1                     deg
  p2                Mechanical power                            W
  torque            Torque                                      Nm
  kt                Torque factor (peak)
  dag               Airgap diameter                             m
  ================  ========================================== =============
  
    Example::

      {'psiq': [0.330062, 0.33031268],
      'psid': [0.08589968, -0.005226678],
      'ld': [0.0008623392, 0.0008623392],
      'lq': [0.0013202480000000002, 0.001458122],
      'psim': [0.08589968, 0.08589968],
      'speed': 50.0,
      'npoles': 8,
      'lfe': 0.08356,
      'psim0': 0.1153128,
      'u1': [145.0, 428.5, 415.1],
      'gamma': [75.44, 90.91],
      'dag': 0.16117,
      'i1': [0, 250.0, 250.0],
      'beta': [0.0, -25.0],
      'kt': [1.14],
      'up0': 145.0,
      'up': 108.0,
      'p2': [80958.54109011423, 127081.80850105156],
      'phi': [50.44, 65.91],
      'torque': [257.69904, 404.51396],
      'cosphi': [0.63688591473536493, 0.40817113454379084]}

Weight
======

  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  total              Total weight                              kg
  conductor          Weight of conductors                      kg
  magnet             Weight of magnets                         kg
  iron               Weight of active iron                     kg
  ================  ========================================== =============

  Example::
    
    {'total': 28.188,
    'iron': 24.165,
    'conductor': 2.853,
    'magnet': 1.17}

Weights
=======

    List of weights (iron, conductors, magnets): in stator and rotor in kg

    Example::
      
    [[18.802, 2.853, 0.0],
    [5.363, 0.0, 1.17],
