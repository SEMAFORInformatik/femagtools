.. _machine-model:

Machine
=======

Three different types of machines are supported:

* RFM: Radial Flux Machines,
* LM: Linear Machines,
* AFPM: Axial Flux Machine.

Machines have a set of basic parameters, a stator, a magnet and a winding:

==============  ======================================  ======
Parameter        Description                            Unit
==============  ======================================  ======
name             Name of machine
lfe              Lenght of iron                         m
afmtype          "S1R1", "S2R1", "S1R2"
poles            Number of poles
outer_diam       Outer diameter (yoke side)             m
bore_diam        Bore diameter  (airgap side)           m
inner_diam       Inner diameter (yoke or shaft)         m
shaft_diam       Shaft diameter                         m
airgap           airgap width                           m
external_rotor   True, False                            False
ffactor          processing factor for iron losses
coord_system     1 (x/y) or 2 (r/z)                     0
dxffile          (see :ref:`model_creation_with_dxf`)
==============  ======================================  ======

.. Note::

   depending on the type (RFM, AFM, LM) not all combinations are useful.
   A LM for example does not have a diameter.

Axial Flux Machines (AFPM)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Axial flux machines have an airgap between stator and rotor that
is aligned parallel with the axis of rotation. Three types of
configurations are distinguished which are defined by the parameter afmtype:

+-----------------------------+-----------------------------+-----------------------------+
| S1R1                        | S2R1                        | S1R2                        |
|                             |                             |                             |
| .. image:: ../afpm/S1R1.svg | .. image:: ../afpm/S1R1.svg | .. image:: ../afpm/S1R2.svg |
|    :width: 200              |    :width: 200              |    :width: 200              |
+-----------------------------+-----------------------------+-----------------------------+

For the simulation the AFPM is split into a number of slices (usually 3) each of which is treated as a linear machine.

.. Note::

   Machine models with any of the above listed types must use the models :ref:`afm_stator`
   and :ref:`afm_rotor`.

.. toctree::
   :maxdepth: 2

   stator
   windings
   magnet
   rotor
   dxf
