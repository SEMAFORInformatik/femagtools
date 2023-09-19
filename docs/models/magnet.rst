
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

   * the mcvkey parameters either reference a filename without extension (Example 'M330-50A') which must be found in the directory defined by the parameter magnetizingCurves of the Femag constructor or the name of an entry in the :ref:`magnetizingCurve` object.
   * the material parameter references a name of the :ref:`magnetMaterial` list.


magnetSector
~~~~~~~~~~~~

  .. image:: ../slot-parameters/magnetSector.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
magn_num                       number of magnets per pole
magn_height     HM      m
magn_width_pct                 magnet width per pole width
magn_width      BM      m      used if magn_width_pct is undefined or 0
condshaft_r             m      conducting shaft radius
magn_rfe                m      inner radius of magnet
magn_len                       relative magnet length
magn_shape              m
bridge_height   BH
bridge_width    BW
magn_ori                       orientation 1: parallel, 2: polar, 3: halbach
magn_type                      1: arc 2: arc par. 3: rect. 4: curved rect.
==============  ======  ====== =============================================

  .. image:: ../slot-parameters/magntype1.png
             :width: 140
  .. image:: ../slot-parameters/magntype2.png
             :width: 140
  .. image:: ../slot-parameters/magntype3.png
             :width: 140
  .. image:: ../slot-parameters/magntype4.png
             :width: 140

magnetIron
~~~~~~~~~~

  .. image:: ../slot-parameters/magnetIron.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
magn_height     HM      m
magn_width      BM      m
gap_ma_iron     DE_FEM  m
air_triangle    BS      m
iron_height     HS      m
magn_rem                T      magn. remanenc
condshaft_r             m      conducting shaft radius if < RI
magn_ori                       orientation 1: parallel, 2: polar, 3: halbach
bridge_height   BH      m
bridge_width    BW      m
iron_shape      HA      m
==============  ======  ====== =============================================

magnetIron3
~~~~~~~~~~~

  .. image:: ../slot-parameters/magnetIron3.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
magn_height     HM      m
magn_width      BM      m
gap_ma_iron     DE_FEM  m
air_triangle    BS      m
iron_height     HS      m
magn_num                       number of magnets
shaft_rad               m      conducting shaft radius if < RI
magn_ori                       orientation 1: parallel, 2: polar, 3: halbach
gap_ma_right    BR      m
gap_ma_left     BL      m
iron_shape      HA      m
==============  ======  ====== =============================================

magnetIron4
~~~~~~~~~~~

  .. image:: ../slot-parameters/magnetIron4.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
magn_height     HM      m
magn_width      BM      m
gap_ma_iron     DE_FEM  m
air_triangle    BS      m
iron_height     HS      m
air_space_h     H_air   m
corner_r        R1      m
magn_dist_ra    DM      m
air_sp_ori                     orientation 0: lin 1: par RA
magn_ori                       orientation 1: parallel, 2: polar, 3: halbach
magn_num                       number of magnets (1 or 2)
iron_shape      HA      m
==============  ======  ====== =============================================

magnetIron5
~~~~~~~~~~~

  .. image:: ../slot-parameters/magnetIron5.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
magn_height     HM      m
magn_width      BM      m
gap_ma_iron     DE_M    m
iron_bfe        BFE     m
iron_height     HS      m
air_space_h     H_air   m
air_space_b     B_air   m
corner_r        R1      m
magn_dist_ra    DM      m
air_sp_ori                     orientation 0: lin 1: par RA
magn_num                       number of magnets (1 or 2)
iron_shape      HA      m
==============  ======  ====== =============================================

magnetIronV
~~~~~~~~~~~

  .. image:: ../slot-parameters/magnetIronV.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
magn_height     HM      m
magn_width      BM      m
magn_angle      ALPHA   Deg
gap_ma_iron     DE_M    m
iron_hs         HS      m
iron_height     BR      m
iron_shape      HA      m
air_triangle    BS      m
condshaft_r             m
magn_num                       number of magnets
==============  ======  ====== =============================================


afm_rotor
~~~~~~~~~

  .. image:: ../slot-parameters/afm_rotor.svg

==============  ======  ====== ============================
Name                    Unit   Comment
==============  ======  ====== ============================
magn_height     hm      m
magn_width      wm             magnet width per pole width
yoke_height     hy      m
==============  ======  ====== ============================


============    ===========================================
Name             Parameter
============    ===========================================
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
it can be used for the model creation as an empty dict (see Note in :ref:`stator_slots_fsl`)::

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


.. include:: userspec.rst
