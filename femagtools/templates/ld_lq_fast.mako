--
-- Ld-Lq- Identification
--
	 
m.move_action     =    0.0 -- rotate 
m.arm_length      =    ${model.get('lfe')*1e3}
m.num_pol_pair    =    m.num_poles/2
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle',0)}
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp')}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.phi_start       =    0.0 --  
m.range_phi       =    720./m.num_poles
m.nu_move_steps   =    ${model.get('num_move_steps')}
m.fc_radius       =    (da1+da2)/4
m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.current         =    ${model['i1_max']}*math.sqrt(2.0)/m.num_par_wdgs
m.num_cur_steps   =    ${model['num_cur_steps']}
m.nu_beta_steps   =    ${model['num_beta_steps']}
m.beta_max        =    ${model['beta_max']}
m.beta_min        =    ${model['beta_min']}

m.pm_eff_aktiv    =    0.0

m.pocfilename    = model..'_'..m.num_poles..'p.poc'

run_models("ld_lq_fast")
