--
-- PM/Rel Losses
--
m.move_action     =    0.0 -- rotate 
m.arm_length      =    ${model.get('lfe')*1e3}
m.num_pol_pair    =    m.num_poles/2
m.speed           =    3000./m.num_pol_pair
m.skew_angle      =    ${model.get('skew_angle',0)}
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp')}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.fc_force_points   =  0.0 --    number move points in air gap       
m.phi_start       =    0.0 --  
m.range_phi       =    720./m.num_poles
m.nu_move_steps   =    ${model.get('num_move_steps')}
m.fc_radius       =    (da1+da2)/4
m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.winding_temp    =    ${model.get('wind_temp')}
m.current         =    1.0
m.ntibfilename    =    '${model.get('ntibfilename')}'

m.pocfilename    = model..'_'..m.num_poles..'p.poc'

run_models("pm_sym_loss")
