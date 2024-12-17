--
-- Torque/Force calculation
--
m.move_action     =    ${model.get('move_action', 0)}
% if model.get('calc_fe_loss', 0):
m.calc_fe_loss    =  ${model['calc_fe_loss']}
% endif
% if  model.get('loss_funct',0):
m.loss_funct      =    ${model.get('loss_funct')}
% endif
% if model.get('lfe',0):
m.arm_length      =    ${model.get('lfe')*1e3}
% endif
% if model.get('move_action', 0) == 0:
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle',0)}
m.fc_radius2      =    0.0
% if model.get('range_start', 0):
m.phi_start       =    ${model['phi_start']}
%endif
% if model.get('range_phi', 0):
m.range_phi       =    ${model['range_phi']}
%endif
% else:
m.speed_linear    =    ${model.get('speed')}
m.skew_linear     =    ${model.get('skew_displ',0)}
m.line            =    0
m.two_pole_wi     =    2*m.pole_width
m.range_x         =    m.two_pole_wi
m.range_y         =    0.0
% if model.get('magn_height', 0):
m.magn_height = ${model['magn_height']*1e3}
% endif
% if model.get('yoke_height', 0):
m.yoke_height = ${model['yoke_height']*1e3}
% endif
if m.yoke_height == nil then
m.yoke_height = 0
end
if m.gap_ma_yoke == nil then
m.gap_ma_yoke = 0
end
m.fc_force_points =  5
m.fcpx_mm1        =   m.npols_gen*m.pole_width +1.0
m.fcpy_mm1        =   -ag/2
m.fcpx_mm2        =   -1.0
m.fcpy_mm2        =   m.fcpy_mm1
m.fcpx_mm3        =   m.fcpx_mm2
if m.model_type == 'S2R1_all' then
m.fcpy_mm3        =   -m.magn_height -m.yoke_height -m.gap_ma_yoke -ag -ag/2
else
m.fcpy_mm3        =   -m.magn_height -m.yoke_height -m.gap_ma_yoke -ag -1
end
m.fcpx_mm4        =   m.fcpx_mm1
m.fcpy_mm4        =   m.fcpy_mm3
m.fcpx_mm5        =   m.fcpx_mm1
m.fcpy_mm5        =   m.fcpy_mm1

% endif
m.nu_force_pat    =    0.0
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp',20.0)}
m.winding_temp    =    ${model.get('wind_temp', 20.0)}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}

m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.current         =    ${model.get('current')}*math.sqrt(2.0)/m.num_par_wdgs
m.angl_i_up       =    ${model.get('angl_i_up', 0)}

m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'
m.period_frac    = ${model.get('period_frac', 1)}

run_models("torq_calc")
