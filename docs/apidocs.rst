Femagtools Package
==================

.. automodule:: femagtools

Basic Classes
-------------
.. autoclass:: Femag
   :members:
      
.. autoclass:: femagtools.model.MachineModel
   :members:

.. autoclass:: femagtools.mcv.MagnetizingCurve
   :members:
      
.. autoclass:: femagtools.fsl.Builder
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
      
