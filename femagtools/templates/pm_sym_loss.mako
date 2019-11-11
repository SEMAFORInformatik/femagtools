--
-- PM/Rel Losses
--
m.move_action     =    0.0 -- rotate
% if 'lfe' in model:
m.arm_length      =    ${model.get('lfe')*1e3}
% endif
m.speed           =    1000.
m.skew_angle      =    ${model.get('skew_angle',0)}
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp')}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.fc_force_points   =  0.0 --    number move points in air gap       
m.phi_start       =    ${model.get('phi_start', 0)}
m.range_phi       =    ${model.get('range_phi', 0)}
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}
m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.winding_temp    =    ${model.get('wind_temp')}
m.current         =    1.0
m.ntibfilename    =    model..'.ntib'
m.period_frac     =    ${model.get('period_frac', 1)}

m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'

run_models("pm_sym_loss")
