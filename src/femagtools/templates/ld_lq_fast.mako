--
-- Ld-Lq- Identification
--
	 
m.move_action     =    0.0 -- rotate
% if isinstance(model, dict) and 'lfe' in model:
m.arm_length      =    ${model.get('lfe')*1e3}
% endif
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle',0)}
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp', 'nil')}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.phi_start       =    ${model.get('phi_start', 0)}
m.range_phi       =    ${model.get('range_phi', 0)}
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}
m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.current         =    ${model['i1_max']}*math.sqrt(2.0)/m.num_par_wdgs
m.num_cur_steps   =    ${model['num_cur_steps']}
m.nu_beta_steps   =    ${model['num_beta_steps']}
m.beta_max        =    ${model['beta_max']}
m.beta_min        =    ${model['beta_min']}

m.pm_eff_aktiv    =    0.0
m.calc_noload     =    ${model.get('calc_noload', 1)}
m.period_frac     =    ${model.get('period_frac', 1)}

m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'
% if model.get('load_ex_cur',0):
m.load_ex_cur    =     ${model.get('load_ex_cur',0)}
run_models("ld_lq_f_cur")
%else:
run_models("ld_lq_fast")
%endif