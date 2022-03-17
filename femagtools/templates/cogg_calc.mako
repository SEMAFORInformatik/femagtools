--
-- Cogging
--
m.move_action     =    0.0 -- rotate
% if model.get('lfe', 0):
m.arm_length      =    ${model.get('lfe')*1e3}
% endif

m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle', 0)}
m.nu_skew_steps   =    ${model.get('num_skew_steps', 0)}
m.magn_temp       =    ${model.get('magn_temp')}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.fc_force_points   =  0.0 --    number move points in air gap       
m.phi_start       =    ${model.get('phi_start', 0)}
m.range_phi       =    ${model.get('range_phi', 0)}
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}

m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}
m.nu_force_pat     =  0.0

m.eval_force1     =    ${model.get('eval_force', 0)}
m.period_frac     =    ${model.get('period_frac', 1)}

m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'
% if model.get('vtu_movie', 0):
m.movie_type = 'vtu'
% else:
m.movie_type = nil
% endif
run_models("cogg_calc")
