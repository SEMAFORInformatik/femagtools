--
-- Cogging
--
% if model.get('lfe', 0):
m.arm_length      =    ${model.get('lfe')*1e3}
% endif
m.move_action     =    ${model.get('move_action', 0)}
% if model.get('move_action', 0) == 0:
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle', 0)}
m.fc_force_points   =  0.0 --    number move points in air gap
m.phi_start       =    ${model.get('phi_start', 0)}
m.range_phi       =    ${model.get('range_phi', 0)}
% else:
m.speed_linear    =    ${model.get('speed')}
m.skew_linear     =    ${model.get('skew_displ',0)}
m.line            =    0
m.two_pole_wi     =    2*m.pole_width
m.range_x         =    m.two_pole_wi
m.range_y         =    0.0
m.fc_force_points =  5
m.fcpx_mm1        =   m.npols_gen*m.pole_width +1.0
m.fcpy_mm1        =   -ag/2
m.fcpx_mm2        =   -1.0
m.fcpy_mm2        =   m.fcpy_mm1
m.fcpx_mm3        =   m.fcpx_mm2
if m.model_type == 'S2R1_all' then
m.fcpy_mm3        =   -m.magn_height -m.yoke_height -m.gap_ma_yoke -ag -ag/2
else
m.fcpy_mm3        =   -m.magn_height -m.yoke_height -m.gap_ma_yoke -1 -ag
end
m.fcpx_mm4        =   m.fcpx_mm1
m.fcpy_mm4        =   m.fcpy_mm3
m.fcpx_mm5        =   m.fcpx_mm1
m.fcpy_mm5        =   m.fcpy_mm1

--m.npols_gen       =  1 -- number of sectors simulated
% endif

m.nu_skew_steps   =    ${model.get('num_skew_steps', 0)}
% if model.get('noload_ex_cur', 0): 
m.nloa_ex_cur = ${model.get('noload_ex_cur',0)}
% endif 
% if model.get("noload_ex_cur", 0): 
m.magn_temp = 20   -- dummy parameter
% else: 
m.magn_temp       =    ${model.get('magn_temp')}
% endif
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
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
