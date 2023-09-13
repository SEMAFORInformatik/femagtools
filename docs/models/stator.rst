**Stator**
----------

Stators have basic parameters and slots:

==============  ===============================  =====================
Parameter        Description                     Default
==============  ===============================  =====================
num_slots        Number of Slots Q
num_slots_gen    Number of Slots in Model        m*Q/gcd(Q, 2p*m)
rlength          Relative iron length            1.0
mcvkey_yoke      Name of lamination material     dummy
mcvkey_teeth     Name of lamination material     dummy
fillfac          stacking factor of lamination   1.0
nodedist         Factor for node distance        1.0
==============  ===============================  =====================

.. Note::

   * if a value for num_slots_gen is missing its value is calculated from the number of slots Q and pole pairs p.
   * the mcvkey parameters either reference a filename without extension (Example 'M330-50A') which must be found in the directory defined by the parameter magnetizingCurves of the Femag constructor or the name of an entry in the :ref:`magnetizingCurve`_ object.

.. _stator:

stator1
~~~~~~~

  .. image:: ../slot-parameters/stator1.svg

==============  ====  ====== =============================================
Name                  Unit   Comment
==============  ====  ====== =============================================
tip_rh1         RH1   m
tip_rh2         RH2   m
slot_rf1        RF1   m
slot_width      SW    m
tooth_width     TW    m
==============  ====  ====== =============================================

==============  ===========================================
Name             Parameter
==============  ===========================================
stator2
                 slot_t1,
                 slot_t2,
                 slot_t3,
                 slot_depth,
                 slot_width,
                 corner_width
==============  ===========================================

.. _statorRotor3:

statorRotor3
~~~~~~~~~~~~

  .. image:: ../slot-parameters/statorRotor3.svg

==============  ====  ====== =============================================
Name                  Unit   Comment
==============  ====  ====== =============================================
slot_height     H     m
slot_h1         H1    m
slot_h2         H2    m
slot_width      SW    m
slot_r1         R1    m
slot_r2         R2    m
wedge_width1    B1    m
wedge_width2    B2    m
middle_line                  0: none, 1: vert. 2: horiz., 3: vert+horiz
tooth_width     TW    m      overwrites B1, B2 if > 0
slot_top_sh                  0: arc 1: line, 2: corner
==============  ====  ====== =============================================

stator4
~~~~~~~

  .. image:: ../slot-parameters/stator4.svg

==============  ====  ====== =============================================
Name                  Unit   Comment
==============  ====  ====== =============================================
slot_height     H     m
slot_h1         H1    m
slot_h2         H2    m
slot_h3         H3    m
slot_h4         H4    m
slot_width      SW    m
slot_r1         R1    m
wedge_width1    B1    m
wedge_width2    B2    m
wedge_width3    B3    m
==============  ====  ====== =============================================

==============  ===========================================
Name             Parameter
==============  ===========================================
statorBG
                 yoke_diam_ins,
                 slot_h1,
                 slot_h3,
                 slot_width,
                 slot_r1,
                 slot_r2,
                 middle_line,
                 tooth_width,
		 tip_rad,
		 slottooth
stator3Linear
                 slot_height,
                 slot_h1,
                 slot_h2,
                 tip_slot,
                 yoke_height,
                 slot_r1,
                 slot_r2,
                 width_bz,
                 tooth_width
afm_stator
                 slot_height,
                 slot_h1,
                 slot_h2,
                 yoke_height,
                 slot_width,
                 slot_open_width,
                 slot_r1,
                 slot_r2
<filename>
                 (see :ref:`stator_slots_fsl`)
dxffile
                 (see :ref:`stator_slots_dxf`)
==============  ===========================================

.. Note::

   All units are metric units.

.. _stator_slots_fsl:

User defined Stator Slots with FSL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a Mako or FSL file that includes the definition of stator geometry exists and is readable it can be used for the model creation.

Example with file mystator.fsl::

  machine = dict(
      name="Motor",
      ...
      stator=dict(
          mcvkey_yoke='dummy',
	  mcvkey_shaft="dummy",
	  mystator=dict()
      ),
      ...

.. Note::
   The file search path can be set with the parameter 'templatedirs' in the Femag or Builder class.

.. _stator_slots_dxf:

User defined Stator Slots with DXF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a DXF file that defines the stator geometry exists and is readable
it can be used to create the FSL of the model.
All DXF conversion parameters are supported.

Example::

  machine = dict(
      name="Motor",
      ...
      stator=dict(
          mcvkey_yoke='dummy',
	  dxffile=dict(
	      name="mystator.dxf",
	      position='out',
              split=True
	  )
      ),
      ...

==========   ============================  =======
Parameters   Description                   Default
==========   ============================  =======
position     'in' or 'out'
split        splits intersecting lines at  False
             their intersection-points
plot         creates the plot              False
	     of the integrated object
==========   ============================  =======

.. Note:: The split option is required only if intersecting lines have no common point.
