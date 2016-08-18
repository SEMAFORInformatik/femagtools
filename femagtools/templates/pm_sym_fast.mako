--
-- PM/Rel Simulation
--
m.move_action     =    0.0 -- rotate 
m.arm_length      =    ${model.get('lfe')*1e3}
m.num_pol_pair    =    m.num_poles/2
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle',0)}
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp')}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.fc_force_points   =  0.0 --    number move points in air gap       
m.phi_start       =    0.0 --  
m.range_phi       =    720./m.num_poles
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}
m.fc_radius       =    (da1+da2)/4
m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.eval_force      =    ${model.get('eval_force', 0)}
m.current         =    ${model.get('current')}*math.sqrt(2.0)/m.num_par_wdgs
m.angl_i_up       =    ${model.get('angl_i_up', 0)}
m.optim_i_up      =    ${model.get('optim_i_up', 0)}

m.pm_eff_aktiv    =    0.0

m.pocfilename    = model..'_'..m.num_poles..'p.poc'

run_models("pm_sym_fast")
