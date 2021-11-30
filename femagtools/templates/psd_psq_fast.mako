--
-- Psid-Psiq- Identification
--
% if model.get('magn_temp',0):
set_dev_data("magn_temp", ${model.get('magn_temp')})
%endif
	 
m.move_action     =    0.0 -- rotate
m.speed           =    ${model.get('speed')*60}
m.skew_angle      =    ${model.get('skew_angle',0)}
m.nu_skew_steps   =    ${model.get('num_skew_steps',0)}
m.fc_mult_move_type =  1.0 --  Type of move path in air gap
m.phi_start       =    ${model.get('phi_start', 0)}
m.range_phi       =    ${model.get('range_phi', 0)}
m.nu_move_steps   =    ${model.get('num_move_steps', 49)}
m.num_par_wdgs    =    ${model.get('num_par_wdgs',1)}

m.maxid           =    ${model['maxid']}/m.num_par_wdgs
m.minid           =    ${model['minid']}/m.num_par_wdgs
m.maxiq           =    ${model['maxiq']}/m.num_par_wdgs
m.miniq           =    ${model['miniq']}/m.num_par_wdgs
m.delta_id        =    ${model['delta_id']}/m.num_par_wdgs
m.delta_iq        =    ${model['delta_iq']}/m.num_par_wdgs
% if model.get('load_ex_cur',0):
m.load_ex_cur     =    ${model['load_ex_cur']}
%endif
m.pm_eff_aktiv    =    0.0
m.calc_noload     =    ${model.get('calc_noload', 1)}
m.period_frac     =    ${model.get('period_frac', 1)}
m.pocfilename    = '${model.get('pocfilename', 'sin.poc')}'
run_models("psd_psq_fast")

