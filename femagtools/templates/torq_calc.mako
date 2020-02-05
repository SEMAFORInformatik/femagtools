--
-- Torque/Force calculation
--
m.move_action     =    ${model.get('move_action', 0)}
% if 'lfe' in model:
m.arm_length      =    ${model.get('lfe')*1e3}
% endif
% if model.get('move_action', 0) == 0:
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle',0)}
m.fc_radius2      =    0.0
m.phi_start       =    ${model.get('phi_start', 0)}
m.range_phi       =    ${model.get('range_phi', 0)}
m.nu_force_pat    =    0.0
% else:
m.speed_linear    =    ${model.get('speed')}
m.skew_linear     =    ${model.get('skew_displ',0)}
m.line            =    0
m.two_pole_wi     =    2*m.pole_width
m.range_x         =    m.two_pole_wi
m.range_y         =    0.0 

m.fc_force_points =  5                                                
m.fcpx_mm1        =   m.npols_gen*m.pole_width +1.0
m.fcpy_mm1        =   -3*ag/4
m.fcpx_mm2        =   -1.0
m.fcpy_mm2        =   m.fcpy_mm1
m.fcpx_mm3        =   m.fcpx_mm2
m.fcpy_mm3        =   -m.magn_height -m.yoke_height -m.gap_ma_yoke -4
m.fcpx_mm4        =   m.fcpx_mm1
m.fcpy_mm4        =   m.fcpy_mm3
m.fcpx_mm5        =   m.fcpx_mm1
m.fcpy_mm5        =   m.fcpy_mm1

m.npols_gen       =  1 -- number of sectors simulated

% endif
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp',20.0)}
m.winding_temp    =    ${model.get('wind_temp', 20.0)}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}

m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.current         =    ${model.get('current')}*math.sqrt(2.0)/m.num_par_wdgs
m.angl_i_up       =    ${model.get('angl_i_up', 0)}

m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'

run_models("torq_calc")
