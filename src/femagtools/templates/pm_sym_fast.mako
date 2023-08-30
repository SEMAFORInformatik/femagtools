--
-- PM/Rel Simulation
--
% if model.get('explicit_mode',0):
set_sim_data("explicit_mode", ${model.get('explicit_mode',0)})
% endif
% if model.get('wind_temp',0):
set_dev_data("cond_temp", ${model.get('wind_temp')}, ${model.get('wind_temp')})
% endif
m.move_action     =    ${model.get('move_action', 0)}
% if  model.get('lfe',0):
m.arm_length      =    ${model.get('lfe')*1e3}
% endif
% if model.get('move_action', 0) == 0:
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle',0)}
m.phi_start       =    ${model.get('phi_start', 0)}
m.range_phi       =    ${model.get('range_phi', 0)}
m.fc_force_points   =  0.0
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

m.npols_gen       =  1 -- nuber of sectors simulated

% endif
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.magn_temp       =    ${model.get('magn_temp')}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}

m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}
m.wdcon           =    ${model.get('wdgcon', 0)}   --    Connection: 0=open, 1=star, 2=delta
m.eval_force      =    ${model.get('eval_force', 0)}
m.current         =    ${model.get('current')}*math.sqrt(2.0)/m.num_par_wdgs
m.angl_i_up       =    ${model.get('angl_i_up', 0)}
m.optim_i_up      =    ${model.get('optim_i_up', 0)}

m.pm_eff_aktiv    =    0.0
m.calc_noload     =    ${model.get('calc_noload', 1)}
m.period_frac     =    ${model.get('period_frac', 1)}

m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'
% if model.get('vtu_movie', 0):
m.movie_type = 'vtu'
% else:
m.movie_type = nil
% endif
run_models("pm_sym_fast")
