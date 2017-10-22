Femagtools Package
==================

.. automodule:: femagtools

Basic Classes
-------------
.. autoclass:: Femag
   :members:
      
.. autoclass:: femagtools.model.MachineModel
   :members:

.. autoclass:: femagtools.fsl.Builder
   :members:

.. automodule:: femagtools.plot
   :members:

Material Handling
-----------------

.. autoclass:: femagtools.mcv.MagnetizingCurve
   :members:

.. autoclass:: femagtools.tks.Reader
   :members:
      
.. autoclass:: femagtools.vbf.Reader
   :members:
      
.. automodule:: femagtools.losscoeffs
   :members:
      

Analytical Models
-----------------
.. autoclass:: femagtools.machine.PmRelMachineLdq
   :members:
      
.. autoclass:: femagtools.machine.PmRelMachinePsidq
   :members:


Calculation Engines
-------------------
.. autoclass:: femagtools.job.Job
   :members:

.. autoclass:: femagtools.job.CloudJob
   :members:
   
.. autoclass:: femagtools.job.Task
   :members:

.. autoclass:: femagtools.job.CloudTask
   :members:

.. autoclass:: femagtools.multiproc.Engine
   :members:

.. automethod:: femagtools.multiproc.run_femag(workdir, fslfile)

.. autoclass:: femagtools.amazon.Engine
   :members:

.. autoclass:: femagtools.google.Engine
   :members:

Parameter Variation and Optimization
------------------------------------
.. autoclass:: femagtools.grid.Grid
.. autoclass:: femagtools.opt.Optimizer
      
