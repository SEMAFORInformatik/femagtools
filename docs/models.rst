Models
******

The models are dictionaries with the properties of a machine or a calculation.

Machine
=======

Machines have a set of basic parameters, a stator, a magnet and a winding:

==============  =============================  ====
Parameter        Description                   Unit
==============  =============================  ====
name             Name of machine
lfe              Lenght of iron                 m
poles            Number of poles
outer_diam       Outer diameter (yoke side)     m
bore_diam        Bore diameter  (airgap side)   m
inner_diam       Inner diameter (yoke)          m
airgap           airgap width                   m
external_rotor   True, False
==============  =============================  ====

Stator
------

Stators have basic parameters and slots:

==============  ============================  =====================
Parameter        Description                  Default  
==============  ============================  =====================
num_slots        Number of Slots Q
num_slots_gen    Number of Slots in Model      m*Q/gcd(Q, 2p*m)
rlength          Relative iron length          1.0
mcvkey_yoke      Name of lamination material
nodedist         Factor for node distance      1.0
==============  ============================  =====================

.. note:: if no value for num_slots_gen is given its value is calculated from the the number of slots Q and pole pairs p.
	  .. version added 0.0.16

Slots
^^^^^
============    ==============  
Name             Parameter      
============    ==============  
stator1  
                 slot_rf1,
                 tip_rh1,
                 tip_rh2, 
                 tooth_width,
                 slot_width
stator2
                 slot_t1,
                 slot_t2,        
                 slot_t3,         
                 slot_depth,      
                 slot_width,      
                 corner_width    
statorRotor3
                 slot_height,
                 slot_h1,    
                 slot_h2,    
                 slot_width, 
                 slot_r1,    
                 slot_r2,
                 wedge_width1,
                 wedge_width2,
                 middle_line, 
                 tooth_width, 
                 slot_top_sh 
statorBG
                 yoke_diam_ins
                 slot_h1,    
                 slot_h3,    
                 slot_width, 
                 slot_r1,    
                 slot_r2,
                 middle_line, 
                 tooth_width,
		 tip_rad,
                 slottooth
============    ==============  

.. note:: all units are metric units

Windings
--------

============    ============================  =======
Name             Parameter                    Default
============    ============================  =======
num_phases      number of phases (m)
num_wires       number of wires per slot
coil_span       coil span
num_layers      number of layers
cufilfact       Fill factor of copper          0.45
culength        rel length of conductor        1.4
slot_indul      insulation thickness in slot   0.0 
============    ============================  =======

Magnet
------

Magnets have basic parameters and slots:

==============  ============================  =======  
Parameter        Description                  Default  
==============  ============================  =======  
mcvkey_yoke      Name of lamination material
mcvkey_shaft     Name of shaft material
material         Name of magnet material
nodedist         Factor for node distance       1.0
==============  ============================  =======

.. note::

  * the mcvkey parameter references a filename without extension (Example 'M330-50A')
  * the material parameter references a name of the `Magnet Material`_ list 

Slots
^^^^^

============    ==============
Name             Parameter      
============    ==============
magnetSector    magn_num,
                magn_width_pct,
                magn_height,
                magn_shape,
                bridge_height,
                magn_type,
                condshaft_r,
                magn_ori,
                magn_rfe,
                bridge_width,
                magn_len
magnetIron      magn_height,
                magn_width,
		gap_ma_iron,
		air_triangle,
		iron_height,
		magn_rem,
		condshaft_r,
		magn_ori,
		bridge_height,
		bridge_width,
		iron_shape
magnetIron2     magn_height,
                magn_width,
		gap_ma_iron,
		air_triangle,
		iron_height,
		magn_rem,
		condshaft_r,
		gap_ma_right,
		gap_ma_left,
		magn_ori,
		iron_shape
magnetIron3     magn_height,
                iron_bfe,
		gap_ma_iron,
		air_triangle,
		iron_height,
		gap_ma_right,
		gap_ma_left,
		condshaft_r,
		magn_num,
		magn_ori,
		iron_shape
magnetIron4     magn_height,
                magn_width,
		gap_ma_iron,
		iron_shape,
		air_space_h,
		iron_bfe,
		magn_di_ra,
		corner_r,
		air_sp_ori,
		magn_ori,
		magn_num
magnetIron5     magn_height,
                magn_width,
		gap_ma_iron,
		iron_bfe,
		air_space_h,
		corner_r,
		air_sp_ori,
		magn_num,
		iron_shape,
		air_space_b,
		magn_di_ra
magnetIronV     magn_height,
                magn_width,
		magn_angle,
		magn_num,
		iron_hs,
		iron_height,
		iron_shape,
		air_triangle,
		gap_ma_iron,
		magn_rem,
		condshaft_r
magnetFC2       yoke_height,
                iron_h1,
		iron_h2,
		iron_b,
		magn_width,
		magn_height,
		iron_bfe,
		iron_bfo,
		iron_shape,
		iron_hp,
		magn_num
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
         mcvkey_shaft="dummy",
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

.. include:: userspec.rst
   
Magnet Material
===============

List of dict objects each having a unique name (or id) and a set of parameters
that describe the magnet properties.

============= ============================   =========   =======
Parameter      Description                    Default     Unit
============= ============================   =========   =======
name          Name of Magnet Material
remanenc      Remanence Induction Br
relperm       Rel. Permeability
spmaweight    Specific Weight                 7500        kg/m³
temcoefbr     Br Temperature Coefficient      -0.001        1/K
temcoefhc     Hc Temperature Coefficient      -0.001        1/K
magntemp      Magnet Temperature              20           °C
magncond      Magnet Conductivity             625000       S m
magnwidth     Magnet Width                    0            m
magnlength    Magnet Length                   0            m
============= ============================   =========   =======


Calculation
===========

Cogging (cogg_calc)

==============  ============================= ==========  ============
Parameter        Description                   Default      Unit
==============  ============================= ==========  ============
speed           Speed                                     1/s
skew_angle      Skewing angle                   0         deg
num_skew_steps  Number of skew steps            0
magn_temp       Magnet Temperature                        °C
num_move_steps  Number of move steps
num_par_wdgs    Number of parallel windings     1      
eval_force      Evaluate force                  0          
==============  ============================= ==========  ============

Example::

  operatingConditions = dict(
    calculationMode="cogg_fast",
    magn_temp=60.0,
    num_move_steps=49,
    speed=50.0)


PM/Rel Machine Simulation (pm_sym_fast)

==============  ============================= ==========  ============
Parameter        Description                   Default      Unit
==============  ============================= ==========  ============
speed           Speed                                     1/s
skew_angle      Skewing angle                   0         deg
num_skew_steps  Number of skew steps            0
magn_temp       Magnet Temperature                        °C
wind_temp       Winding Temperature             20        °C
num_move_steps  Number of move steps
num_par_wdgs    Number of parallel windings     1      
eval_force      Evaluate force                  0         
current         Phase current                             A (RMS)
angl_i_up       Angle I vs. Up                  0         deg
optim_i_up      Optimize Current                0
==============  ============================= ==========  ============

Example::

  operatingConditions = dict(
    calculationMode="pm_sym_fast",
    wind_temp=60.0,
    magn_temp=60.0,
    current=50.0,
    speed=50.0)
  
Ld-Lq Identification (ld_lq_fast)

==============  ============================= ==========  ============
Parameter        Description                   Default      Unit
==============  ============================= ==========  ============
speed           Speed                                     1/s
skew_angle      Skewing angle                   0         deg
num_skew_steps  Number of skew steps            0
magn_temp       Magnet Temperature                        °C
num_move_steps  Number of move steps
num_par_wdgs    Number of parallel windings     1      
eval_force      Evaluate force                  0         
i1_max          Max. phase current                        A (RMS)
beta_min        Min. Beta angle                           deg
beta_max        Max. beta angle                           deg
num_cur_steps   Number of current steps
num_beta_steps  Number of beta steps
==============  ============================= ==========  ============

Example::

  feapars = dict(
    num_move_steps=25,
    calculationMode="ld_lq_fast",
    magn_temp=60.0,
    i1_max=150.0,
    beta_max=0.0,
    beta_min=-60.0,
    num_cur_steps=3,
    num_beta_steps"=4,
    speed=50.0)
  

Psid-Psiq Identification (psd_psq_fast)

==============  ============================= ==========  ============
Parameter        Description                   Default      Unit
==============  ============================= ==========  ============
speed           Speed                                     1/s
skew_angle      Skewing angle                   0         deg
num_skew_steps  Number of skew steps            0
magn_temp       Magnet Temperature                        °C
num_move_steps  Number of move steps
num_par_wdgs    Number of parallel windings     1      
eval_force      Evaluate force                  0         
max_id          Max. Amplitude Id current                 A 
min_id          Min. Amplitude Id current                 A 
max_iq          Max. Amplitude Iq current                 A 
min_iq          Min. Amplitude Iq current                 A 
delta_id        Delta of Id current steps                 A
delta_iq        Delta of Iq current steps                 A
==============  ============================= ==========  ============

Example::

  feapars = dict(
    num_move_steps=25,
    calculationMode="psd_psq_fast",
    magn_temp=60.0,
    max_id=0.0,
    min_id=-150.0,
    max_iq=150.0
    min_iq=0.0,
    delta_id=50.0,
    delta_iq=50.0,
    speed=50.0)
