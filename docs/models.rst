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

Stators have basic parameters and slots:

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

Windings
--------

============    ==========================================
Name             Parameter      
============    ==========================================
num_phases      number of phases
num_wires       number of wires per slot
coil_span       coil span
num_layers      number of layers
cufilfact       Fill factor of copper (default=0.45)
culength        rel length of conductor (default=1.4)
slot_indul      insulation thickness in slot (default=0.0) 
============    ==========================================

Magnet
------

Magnets have basic parameters and slots:

==============  ============================  ====
Parameter        Description                  Unit
==============  ============================  ====
mcvkey_yoke      Name of lamination material
mcvkey_mshaft    Name of shaft material
==============  ============================  ====

Slots
^^^^^

============    ==============
Name             Parameter      
============    ==============
magnetSector    magn_num
                magn_width_pct
                magn_height
                magn_shape
                bridge_height
                magn_type
                condshaft_r
                magn_ori
                magn_rfe
                bridge_width
                magn_len
============    ==============

Example::
  
  machine = dict(
     name="PM 130 L4",
     lfe=0.1,
     poles=4,
     outer_diam=0.13,
     bore_diam=0.07,
     inner_diam=0.015,
     airgap=0.001,
     
     stator=dict(
         num_slots=12,
         num_slots_gen=3,
         mcvkey_yoke="dummy",
         rlength=1.0,
         stator1=dict(
             slot_rf1=0.057,
             tip_rh1=0.037,
             tip_rh2=0.037,
             tooth_width=0.009,
             slot_width=0.003)
	 ),

     magnet=dict(
         mcvkey_mshaft="dummy",
         mcvkey_yoke="dummy",
         magnetSector=dict (
	     magn_num=1,
	     magn_width_pct=0.8,
	     magn_height=0.004,
	     magn_shape=0.0,
	     bridge_height=0.0,
	     magn_type=1,
	     condshaft_r=0.02,
	     magn_ori=2,
	     magn_rfe=0.0,
	     bridge_width=0.0,
	     magn_len=1.0 )
	 ),

      windings=dict(
           num_phases=3,
           num_wires=100,
           coil_span=3.0,
           num_layers=1)
  )
  
 

Calculation
===========

Cogging (cogg_calc)

==============  =======================================  ============
Parameter        Description                             Unit
==============  =======================================  ============
speed           Speed                                    1/s
skew_angle      Skewing (default 0)                      deg
num_skew_steps  Number of skew steps (default 0)
magn_temp       Magnet Temperature                       deg Celsius
num_move_steps  Number of move steps
num_par_wdgs    Number of parallel windings (default 1)
eval_force      Evaluate force  (default False)          True, False
==============  =======================================  ============

Example::

  operatingConditions = dict(
    calculationMode="cogg_fast",
    magn_temp=60.0,
    num_move_steps=49,
    speed=50.0)


PM/Rel Machine Simulation (pm_sym_fast)

Ld-Lq Identification (ld_lq_fast)

Psid-Psiq Identification (psd_psq_fast)
