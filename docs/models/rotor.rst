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

   * the mcvkey parameters either reference a filename without extension (Example 'M330-50A') which must be found in the directory defined by the parameter magnetizingCurves of the Femag constructor or the name of an entry in the :ref:`magnetizingCurve` object.

rot_hsm
~~~~~~~

  .. image:: ../slot-parameters/rot_hsm.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
gap_pol_shaft           m
core_height     HPK     m
pole_height     HPO     m
pole_rad        RPS     m
core_width2     BPK2    m
core_width1     BPK     m
pole_width_r    BPSS    m
pole_width      BPSU    m
slot_width      BDNS    m
slot_height     HDNS    m
damper_diam     DD      m
damper_div      TAUD    m
==============  ======  ====== =============================================

rotorAsyn
~~~~~~~~~

  .. image:: ../slot-parameters/asynRotor.svg

==============  ======  ====== =============================================
Name                    Unit   Comment
==============  ======  ====== =============================================
slot_bs2        BS2     m
slot_hs2        HS2     m
slot_b32        B32     m
slot_h32        H32     m
slot_b42        B42     m
slot_h42        H42     m
slot_b52        B52     m
slot_b62        B62     m
slot_h52        H52     m
slot_h62        H62     m
slot_h72        H72     m
==============  ======  ====== =============================================

============  ===========================================
Name             Parameter
============  ===========================================
statorRotor3  (same as in :ref:`statorRotor3`)
rotorKs2      slot_angle, slot_height, slot_topwidth.
              slot_width, slot_h1,
              slot_h2, slot_r1, slot_r2, middle_line
============  ===========================================
