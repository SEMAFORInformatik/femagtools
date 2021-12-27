--
-- SM Simulation with excitation current
--
m.wdcon           = ${model.get('wdgcon', 1)}   --    Connection: 0=open, 1=star, 2=delta
m.move_action     = ${model.get('move_action', 0)}
m.speed           = ${model.get('speed')*60}
m.skew_angle      = ${model.get('skew_angle',0)}
m.nu_skew_steps   = ${model.get('num_skew_steps',0)}
m.eval_force      = ${model.get('eval_force', 0)}
m.num_par_wdgs    = ${model.get('num_par_wdgs',1)}
m.current         = ${model.get('current')}*math.sqrt(2.0)/m.num_par_wdgs
m.angl_i_up       = ${model.get('angl_i_up', 0)}

m.optim_i_up      = ${model.get('optim_i_up', 0)}
m.nu_move_steps   = ${model.get('num_move_steps', 49)}
m.range_phi       = ${model.get('range_phi', 0)}
m.phi_start       = ${model.get('phi_start', 0)}
m.pm_eff_aktiv    = 0   --    Generate additional output (>=1)
m.fc_mult_move_type = 1  --    Type of move path in air gap (1=circ)
m.calc_noload     =    ${model.get('calc_noload', 1)}
m.period_frac    = ${model.get('period_frac', 1)}
m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'
% if model.get('vtu_movie', 0):
m.movie_type = 'vtu'
%endif
-- Excitation current
m.nloa_ex_cur   = ${model.get('nload_ex_cur', 0)} -- No Load Exciting current
m.load_ex_cur   = ${model.get('load_ex_cur', 0)} -- Load Exciting current
m.wdgkeyex      =  0 -- automatic winding selection

run_models("pm_sym_f_cur")
