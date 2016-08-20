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
weight             Total weight and weight of iron, conductor and magnets
================  =======================================================

* Flux: list of dictionaries for each winding

  ================  =======================================================
  Attribute          Description     
  ================  =======================================================
  displ             position
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
  a                 amplitude of sin 
  b                 amplitide of cos 
  ================  =======================================================

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
  ================  =======================================================
