
-- M-Sector-Lin

m.l_corner_x0     = 0.0
m.l_corner_y0     = 0.0
m.magn_height     = ${model['magn_height']*1e3}
% if model.get('magn_width_pct'):
m.magn_width      = ${model['magn_width_pct']*1e2}
% else:
m.magn_width      = -${model['magn_width']*1e3}
% endif
m.pole_width      = ${model['pole_width']*1e3}
m.yoke_heigth     = ${model['yoke_height']*1e3}
m.magn_len        = ${model['magn_len']*100}
m.magn_rem        = m.remanenc
m.gap_ma_yoke     = ${model['gap_ma_yoke']*1e3}
m.magn_ori        = ${model['magn_ori']}
m.airgap_shape    = ${model['airgap_shape']*1e3}
m.magn_type       = ${model['magn_type']}

m.zeroangl        = ${model.get('zeroangle',0)}
 
 m.mcvkey_yoke =   mcvkey_yoke
 
 pre_models("M-Sector-Lin");
