
-- Commutator-Motor-Simulation             
m.curr_step  =  ${model.get('curr_step',1)} --Number of current steps
m.current    =  ${model.get('current',1)} --Nominal winding current [peak] IN [A]
m.brush_angl =  ${model.get('brush_angl',0)}
m.skew_angle = ${model.get('skew_angle',0)}
m.num_skew_st  = ${model.get('num_skew_st',0)} --Number of skew sections: 0=infinite
m.brush_width  = ${model.get('brush_width',1e-2)} -- Rel. brush width (Ref:taupe:piD/2

m.num_arm_bar  = ${model.get('num_arm_bar',0)}  --ZA : number of armature bars
m.num_fi_turn  = ${model.get('num_fi_turn',0)}  --ZF : number of field windings
m.num_par_wdgs = ${model.get('num_par_wdgs',2)}  --N.par.arm.circ: LapW:2a=2p WaveW:2a=2
m.ex_condition = ${model.get('ex_condition',1)}
m.speed        = ${model.get('speed',25)*60}
m.magn_temp    = ${model.get('magn_temp',20)}
m.calc_noload  = ${model.get('calc_noload',1)}
m.forcedens = ${model.get('forcedens',0)}
m.calc_react = ${model.get('calc_react',0)}
m.calc_fe_los = ${model.get('calc_fe_los',0)}
m.delta_br_an = ${model.get('delta_br_an',0)}
m.nu_move_steps = ${model.get('num_move_steps',49)}
m.cmmfilename = model .. "_" .. m.num_poles .."p.cmm"
m.range_phi = ${model.get('range_phi',0)}
m.phi_start = 0

run_models("com-motor_sim")
