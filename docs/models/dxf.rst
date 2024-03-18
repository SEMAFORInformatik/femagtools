
.. _model_creation_with_dxf:

Model Creation with DXF
=======================

The goal of the dxfsl modules is to create a complete FE model for rotating PM and Reluctance Machines on the basis of a DXF file with
as little restrictions as possible.
This has the consequence that symmetries, subregions, magnets, windings, boundary conditions need to be identified.

The procedure is as follows:

1. Read the DXF file and create a graph object with nodes and edges.
2. Identify the areas stator, rotor and airgap and their symmetry axis.
3. Add auxiliary lines if required by the meshing process.
4. Convert the graph object into FSL code including the identified subregions.

For monitoring and trouble-shooting purposes it is possible to create plots that display the intermediate results.

Example with a single dxf file motor.dxf::

   machine = dict(
      name="Motor",
      lfe=0.001,

      dxffile=dict(
         name='motor.dxf'
      ),
      stator=dict(
         mcvkey_yoke='dummy',
	 mcvkey_shaft="dummy"
      ),
      magnet=dict(
         mcvkey_yoke="dummy",
	 mcvkey_shaft="dummy"
      ),
      ...

.. Note::
    * The parameters *poles*, *outer_diam*, *bore_diam* and *airgap* as well as *num_slots* and *num_slots_gen* will be set automatically
    * for additional keys of dxffile see :ref:`stator_slots_dxf`.

Example with two separate dxf files for stator and rotor::

   machine = dict(
      name="Motor",
      lfe=0.001,

      stator=dict(
          mcvkey_yoke='dummy',
  	  dxffile=dict(
	      name='mystator.dxf',
	      position='out'
	 )
      ),
      magnet=dict(
          mcvkey_yoke="dummy",
	  mcvkey_shaft="dummy",
	  dxffile=dict(
	      name='myrotor.dxf',
              position='in'
	  )
      ),
      ...
