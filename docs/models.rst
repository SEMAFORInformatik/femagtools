Models
******

The models are dictionaries with the properties of a machine or a calculation.

Machine
=======

Machines have a set of basic parameters, a stator, a magnet and a winding:

==============  =================  ====
Parameter        Description       Unit
==============  =================  ====
name             Name of machine
lfe              Lenght of iron     m
poles            Number of poles
outer_diam       Outer diameter     m
bore_diam        Bore diameter      m
inner_diam       Inner diameter     m
airgap           airgap width       m
external_rotor   True, False
==============  =================  ====

Stator
------

Stators basic parameters and slots:

==============  ============================  ====
Parameter        Description                  Unit
==============  ============================  ====
num_slots        Number of Slots
num_slots_gen    Number of Slots in Model
rlength          Relative iron length
mcvkey_yoke      Name of lamination material
==============  ============================  ====


Slots
^^^^^
============    ==============  
Name             Parameter      
============    ==============  
stator1  
                 slot_rf1
                 tip_rh1 
                 tip_rh2 
                 tooth_width
                 slot_width
stator2
                 slot_t1
                 slot_t2         
                 slot_t3         
                 slot_depth      
                 slot_width      
                 corner_width    
statorRotor3
                 slot_height
                 slot_h1    
                 slot_h2    
                 slot_width 
                 slot_r1    
                 slot_r2
                 wedge_width1
                 wedge_width2
                 middle_line 
                 tooth_width 
                 slot_top_sh 
============    ==============  

Calculation
===========

Cogging (cogg_calc)

==============  ============================  ============
Parameter        Description                  Unit
==============  ============================  ============
speed           Speed                         1/s
skew_angle      Skewing                       deg
num_skew_steps  Number of skew steps
magn_temp       Magnet Temperature            deg Celsius
num_move_steps  Number of move steps
num_par_wdgs    Number of parallel windings
eval_force1     Evaluate force                True, False
==============  ============================  ============

PM/Rel Machine Simulation (pm_sym_fast)

Ld-Lq Identification (ld_lq_fast)

Psid-Psiq Identification (psd_psq_fast)
