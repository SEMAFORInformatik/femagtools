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
windings           Winding properties
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
inertia            Inertia of stator and rotor
leak_dist_wind     End Winding Leakage
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

* torque: list of dictionaries for each current and/or beta angle (No load, load current with beta=0 and load current with beta)

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  angle             Position
  current_1         Current
  force_x           Force in x-direction 
  force_y           Force in y-direction 
  t_idpsi           Torque with dq-parameters
  torque            Torque with Maxwell stress tensor
  ripple            Diff between max and min torque value
  ================  =======================================================

  .. Note:: the force values are valid for the simulated model segment only.
	    The torque values are valid for the complete machine.

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

  ================  =============================  ========================
  Attribute          Description                   Unit
  ================  =============================  ========================
  iq                Iq current list (n)            A
  id                Id current list (m)            A
  psid              Psid matrix (n x m)            Vs
  psiq              Psiq matrix (n x m)            Vs
  torque            Torque matrix (n x m)          Nm
  losses            dict of loss values          
  ================  =============================  ========================

  * losses
    
  ================  ====================================  =====
  Attribute          Description                          Unit
  ================  ====================================  =====
  styoke            Losses of stator yoke (n x m)         W 
  stteeth           Losses of stator teeth (n x m)        W
  rotor             Losses of rotor (n x m)               W
  magnet            Losses of magnet (n x m)              W
  styoke_hyst       Hyst. Losses of stator yoke (n x m)   W
  styoke_eddy       Eddy Losses of stator yoke (n x m)    W
  stteeth_hyst      Hyst. Losses of stator teeth (n x m)  W
  stteeth_eddy      Eddy Losses of stator yoke (n x m)    W
  rotor_hyst        Hyst. Losses of rotor (n x m)         W
  rotor_eddy        Eddy Losses of rotor (n x m)          W
  speed             Speed                                 1/s
  ================  ====================================  =====
    
Psidq Ldq
=========

  ================  =============================  ========================
  Attribute          Description                   Unit
  ================  =============================  ========================
  iq                Iq current list (n)            A
  id                Id current list (m)            A
  ld                Ld matrix (n x m)              H
  lq                Lq matrix (n x m)              H
  psim              Psim matrix (n x m)            Vs
  psid              Psid matrix (n x m)            Vs
  psiq              Psiq matrix (n x m)            Vs
  torque            Torque matrix (n x m)          Nm
  ================  =============================  ========================

Ldq
===

  ================  =============================  ========================
  Attribute          Description                   Unit
  ================  =============================  ========================
  i1                I1 current list (n)            A
  beta              Beta current angle list (m)    deg
  ld                Ld matrix (n x m)              H
  lq                Lq matrix (n x m)              H
  psim              Psim matrix (n x m)            Vs
  psid              Psid matrix (n x m)            Vs
  psiq              Psiq matrix (n x m)            Vs
  torque            Torque matrix (n x m)          Nm
  ================  =============================  ========================

  
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
  p_sim             Number of poles in model
  Q                 Total number of stator slots
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
  pocfile           Name of POC file used                 
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
  u1_sim            Terminal voltage  (Sim)                     V
  u1_fe             Terminal voltage  (FE)                      V
  up                MMF voltage                                 V
  up0               MMF voltage at-noload                       V
  gamma             Angle between Up and U1                     deg
  i1                Phase current                               A
  phi               Angle between U1 and I1                     deg
  p2                Mechanical power                            W
  torque            Torque                                      Nm
  torque_sim        Torque (Sim)                                Nm
  torque_fe         Torque (FE)                                 Nm
  kt                Torque factor (peak)
  dag               Airgap diameter                             m
  ================  ========================================== =============
  
    Example::

      {'i1': [0, 243.3, 243.3],
      'beta': [0.0, -35.54],
      'ld': [0.0005299380000000001, 0.0005299380000000001],
      'lq': [0.0012425400000000003, 0.0014455000000000002],
      'torque': [444.97800000000007, 829.0680000000001],
      'kt': [2.41],
      'psim0': 0.10266,
      'up0': 258.0,
      'psim': [0.10320280000000001, 0.10320280000000001],
      'speed': 66.66666666666667,
      'npoles': 12,
      'lfe': 0.11800000000000001,
      'dag': 0.3132,
      'u1': [258.0, 805.1564407729971, 727.4119501308116],
      'gamma': [70.67427472336531, 83.94107114196522],
      'phi': [70.67427472336531, 48.401071141965225],
      'cosphi': [0.33093811730811373, 0.6639122324852806],
      'psid': [0.10320280000000001, 0.028249200000000002],
      'psiq': [0.30231600000000003, 0.28615],
      'torque_fe': [452.0, 836.0],
      'torque_sim': [444.9, 829.0],
      'p2': [186391.9487745439, 347279.1917501844],
      'u1_fe': [801.6, 722.6],
      'u1_sim': [802.9, 722.6],
      'up': [259.3769266479174, 259.3769266479174]}

Magnet
======

  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  Br                 Remanence                                 T 
  Hc                 Coercitivity                              A/m
  muer               rel Permeability                            
  Tmag               Temperature                               °C
  alpha              Temperature coefficient of Br             1/K
  demag_pc           Demag Limit                               %
  demag_hx           Demag Limit                               A/m
  area               Area                                      mm²
  sigma_PM           El. Conductivity                          S/m
  ================  ========================================== =============

  Example::
    
    {Br': 1.2,
    'Hc': -909.457,
    'muer': 1.05,
    'Tmag': 120.0,
    'alpha': -0.1,
    'demag_pc': 95.0,
    'demag_hx': -863.984,
    'area': 4136.087,
    'sigma_PM': 625000.0}

    
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

Inertia
=======

    List of inertia (Stator, rotor) [Unit kg m²/mm]

    Example::
      
       [0.23, 0.39]

Windings
========

  Dictionary with winding key:
  
  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  dir                list of winding directions 
  N                  list with number of conductors
  R                  list of radius                            m
  PHI                list of angles                            deg
  ================  ========================================== =============

  Example::

     {  1: {'N': [4.0, 4.0, 4.0, 4.0],
            'R': [92e-3, 92.0086, 92e-3, 92e-3],
            'dir': [1, 1, 1, -1],
            'PHI': [3.0203, 4.4797, 11.9797, 40.5202]},
        2: {'N': [4.0, 4.0, 4.0, 4.0],
            'R': [92e-3, 92e-3, 92e-3, 92.0086],
            'dir': [1, 1, 1, 1],
            'PHI': [25.5202, 33.0202, 34.4797, 41.9797]},
        3: {'N': [4.0, 4.0, 4.0, 4.0],
            'R': [92e-3, 92e-3, 92e-3, 92e-3],
            'dir': [-1, -1, -1, -1],
            'PHI': [10.5202, 18.0202, 19.4797, 26.9797]}
     }
 
Losses
======

 List of dictionaries with losses for noload and load calculation:

  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  beta               angle I Up                                °
  current            winding current (RMS)                     A
  staza              losses in stator teeth                    W
  stajo              losses in stator yoke                     W
  rotfe              losses in rotor                           W
  winding            losses in windings                        W
  magnetB            losses in magnet (B-Method)               W
  magnetJ            losses in magnet (J-Method)               W
  total              total losses                              W
  r1                 winding resistance                        Ohm
  fft                dict of harmonic spectrum losses
                     rotor, staza, stajo
		     with: order, freq, hyst, eddy
  ================  ========================================== =============

  Example::
    
    {'beta': 0.0,
     'current': 0.0,
     'magnetB': 0.0,
     'magnetJ': 0.053,
     'r1': 0.0,
     'rotfe': 483.806,
     'stajo': 1242.913,
     'staza': 1664.52,
     'total': 3391.292,
     'winding': 0.0,
     'fft': {
        'rotor': {'eddy': (475.937,),
                  'freq': (600.0,),
                  'hyst': (7.869,),
                  'order_el': (6,)},
        'stajo': {'eddy': (15.777, 138.777, 394.427, 206.489, 313.927, 134.139, 5.329),
                  'freq': (100.0, 300.0, 500.0, 700.0, 900.0, 1100.0, 1300.0),
                  'hyst': (9.983, 9.178, 9.391, 2.508, 2.307, 0.66, 0.019),
                  'order_el': (1, 3, 5, 7, 9, 11, 13)},
        'staza': {'eddy': (13.06, 117.544, 325.934, 212.208, 417.321, 326.528, 220.231),
                  'freq': (100.0, 300.0, 500.0, 700.0, 900.0, 1100.0, 1300.0),
                  'hyst': (8.135, 7.774, 7.76, 2.578, 3.067, 1.606, 0.776),
                  'order_el': (1, 3, 5, 7, 9, 11, 13)}}
     }


Demag
=====

 List of dictionaries with demagnetization information

  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  displ              rotor position                            °
  current            winding current (RMS)                     A
  beta               angle I Up                                °
  current_1          current winding 1 (RMS)                   A
  current_2          current winding 2 (RMS)                   A
  current_3          current winding 3 (RMS)                   A
  H_max              maximum field strength                    kA/m
  H_av               average field strength                    kA/m
  area               area with H > Hx                          %
  ================  ========================================== =============

Leak_dist_wind
==============

  Dict with end-winding leakage values (version added 0.9.9)

  ================  ========================================== =============
  Attribute          Description                               Unit
  ================  ========================================== =============
  nseg              Number of segments
  npolsim           Number of poles in model
  fc_radius         Force radius (center of airgap             m
  armatureLength    Lenght of armature                         m
  perimrad          Radius of perimeter                        m
  vbendrad          Bending radius vertical                    m
  endheight         End winding height                         m
  wiredia           Diameter of wire                           m
  L0e               Ext. Inductance                            H
  Lde               Ext. Inductance                            H
  Lqe               Ext. Inductance                            H
  L0i               Int. Inductance                            H
  Ldi               Int. Inductance                            H
  Lqi               Int. Inductance                            H
  ================  ========================================== =============

scData (Short Circuit)
======================

  Dict  values of short circuit calculation (version added 0.9.30)

  ===================    =================================== =========
  Attribute              Description                            Unit
  ===================    =================================== =========
  speed                  Speed                                   1/s
  ikd                    stationary phase current amplitude      A
  tkd                    stationary Torque                       Nm
  iks                    Peak Current                            A
  tks                    Peak Torque                             Nm
  time                   Time vector                             s
  ia                     Phase a Current vector                  A
  ib                     Phase b Current vector                  A
  ic                     Phase c Current vector                  A
  peakWindingCurrents    peak current of each phase              A
  ===================    =================================== =========
