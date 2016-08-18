Introduction and Overview
*************************

Run FEMAG with FSL
++++++++++++++++++
Example (single process)::
  
  workdir = os.path.join(os.environ['HOME'], 'femag')
  femag = femagtools.Femag(workdir)
  femag.run('femag.fsl')

multi processes::

  engine = femagtools.MultiProc()
  job = engine.create_job(workdir)
  for fsl in ('femag-1.fsl', 'femag-2.fsl', 'femag-3.fsl'):
      task = job.add_task()
      task.add_file(fsl)

  numtasks = engine.submit()
  status = engine.join()
  
Read BCH/BATCH File
+++++++++++++++++++
Example::

  bch = femagtools.read_bchfile('TEST_002.BCH')
  print(bch.machine['torque'])


Create FSL from Templates
+++++++++++++++++++++++++
Example::

  machine = dict(
     name = "PM 130 L4",
     lfe = 0.1,
     poles = 4,
     outer_diam = 0.13,
     bore_diam = 0.07,
     inner_diam = 0.015,
     airgap = 0.001,
     
     stator = dict(
         num_slots = 12,
         num_slots_gen = 3,
         mcvkey_yoke = "dummy",
         rlength = 1.0,
         stator1 = dict(
             slot_rf1 = 0.057,
             tip_rh1 = 0.037,
             tip_rh2 = 0.037,
             tooth_width = 0.009,
             slot_width = 0.003)
	 ),

     magnet = dict(
         mcvkey_mshaft = "dummy",
         mcvkey_yoke = "dummy",
         magnetSector = dict (
	     magn_num = 1,
	     magn_width_pct = 0.8,
	     magn_height = 0.004,
	     magn_shape = 0.0,
	     bridge_height = 0.0,
	     magn_type = 1,
	     condshaft_r = 0.02,
	     magn_ori = 2,
	     magn_rfe = 0.0,
	     bridge_width = 0.0,
	     magn_len = 1.0 )
	 ),

      windings = dict(
           num_phases = 3,
           num_wires = 100,
           coil_span = 3.0,
           num_layers = 1)
  )
  builder = femagtools.FslBuilder()
  model = femagtools.Model(machine)
  fsl = builder.create_model(model)
  with open('femag.fsl', 'w') as f:
      f.write('\n'.join(fsl))

After opening this file in FEMAG the shown geometry is created:

.. image:: geom.png
   :height: 240pt

PM machine characteristics
++++++++++++++++++++++++++

Definition of PM machine with Ld,Lq parameters::

  p = 4
  r1 = 0.0806
  le = 0.0
  ls = 0.0
  wind_temp = 20.0
  ld = [0.0014522728, 0.0014522728]
  lq = [0.0032154, 0.0038278836]
  psim = [0.11171972000000001, 0.11171972000000001]
  i1 = [80.0]
  beta = [0.0, -41.1]

  pm = femagtools.PmRelMachineLdq(3, p,
                                  psim,
                                  ld,
                                  lq,
                                  r1,
                                  beta,
                                  i1)

Calculation of minimal current and frequency at given torque and max voltage::

  tq = 170.0
  u1 = 340.0

  iqx, idx = pm.iqd_torque(tq)
  w1 = pm.w1_u(u1, idx, iqx)
  i1 = np.linalg.norm(np.array((iqx, idx)))

.. plot:: pyplots/pmfieldweak.py
      
Speed-Torque characteristics with max power::

  def torque(T, pmax, wm):
      """shaft torque as a function of rotor angular speed"""
      if wm <= pmax / T:
          return T
      return pmax / wm


  pmax = 60e3
  n = np.linspace(0, 75, 20)
  T = [torque(Tmax, pmax, 2*np.pi*nx) for nx in n]
  r = pm.characteristics(T, n, u1)

.. plot:: pyplots/pmchar.py
  

Multi-Objective Optimization
++++++++++++++++++++++++++++

Example::

  engine = femagtools.Condor()
  opt = femagtools.Optimizer(parameters, engine, workdir)
  num_generations = 10
  results = opt.optimize(num_generations)
