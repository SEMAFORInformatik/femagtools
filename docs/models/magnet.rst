
**Magnet**
----------

Magnets have basic parameters and slots:

==============  ================================  =======
Parameter        Description                      Default
==============  ================================  =======
mcvkey_yoke      Name of lamination material      dummy
mcvkey_shaft     Name of shaft material           dummy
material         Name of magnet material
nodedist         Factor for node distance         1.0
==============  ================================  =======

.. Note::

   * the mcvkey parameters either reference a filename without extension (Example 'M330-50A') which must be found in the directory defined by the parameter magnetizingCurves of the Femag constructor or the name of an entry in the `magnetizingCurve`_ object.
   * the material parameter references a name of the `Magnet Material`_ list.

Magnet Slots
^^^^^^^^^^^^

============    ===========================================
Name             Parameter
============    ===========================================
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
afm_rotor       yoke_height
                magn_height
                magn_width
                spoke_width

<filename>      see :ref:`rotor_slots_fsl`
dxffile         see :ref:`rotor_slots_dxf`
============    ===========================================

**Example**::

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

.. _rotor_slots_fsl:

User defined Magnet Slots with FSL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Example**

If a Mako or FSL file that creates the magnet geometry exists and is readable
it can be used for the model creation as an empty dict (see Note in `stator_slots_fsl`_)::

  machine = dict(
      name="Motor",
      ...
      magnet=dict(
          mcvkey_yoke='dummy',
	  mcvkey_shaft="dummy",
	  myrotor=dict()
      ),
      ...

.. _rotor_slots_dxf:

User defined Magnet Slots with DXF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a DXF file that defines the magnet geometry exists and is readable
it can be used to create the FSL for the model.

Example::

  machine = dict(
      name="Motor",
      ...
      magnet=dict(
          mcvkey_yoke='dummy',
	  mcvkey_shaft="dummy",
	  dxffile=dict(
	      name='mymagnet.dxf',
	      position='in',
              split=True
	  )
      ),
      ...

==========   ============================  =======
Parameters   Description                   Default
==========   ============================  =======
position     'in' or 'out'
split        splits intersecting lines at  False
             their intersection points
plot         creates the plot              False
	     of the integrated object
==========   ============================  =======

.. Note:: The split option is required only if intersecting lines have no common point.

**Rotor**
---------
Rotors have a excitation or short circuit cage winding
and following basic parameters and slots:

==============  ================================  =======
Parameter        Description                      Default
==============  ================================  =======
mcvkey_yoke      Name of lamination material      dummy
mcvkey_shaft     Name of shaft material           dummy
fillfac          stacking factor of lamination    1.0
nodedist         Factor for node distance         1.0
==============  ================================  =======

.. Note::

   * the mcvkey parameters either reference a filename without extension (Example 'M330-50A') which must be found in the directory defined by the parameter magnetizingCurves of the Femag constructor or the name of an entry in the `magnetizingCurve`_ object.

.. include:: userspec.rst
