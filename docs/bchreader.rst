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
