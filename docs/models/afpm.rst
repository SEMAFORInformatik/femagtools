
.. _afpm:
   
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
