
-- M-Sector-Lin

m.l_corner_x0     = 0.0
m.l_corner_y0     = 0.0
m.magn_height     = ${model.get(['magnet', 'magnetSectorLinear', 'magn_height'])*1e3}
<<<<<<< HEAD
m.magn_width      = ${model.get(['magnet', 'magnetSectorLinear', 'magn_width_pct'])*100}
=======
% if model.get(['magnet', 'magnetSectorLinear', 'magn_width_pct']):
m.magn_width      = ${model.get(['magnet', 'magnetSectorLinear', 'magn_width_pct'])*1e2}
% else:
m.magn_width      = -${model.get(['magnet', 'magnetSectorLinear', 'magn_width'])*1e3}
% endif
>>>>>>> development
m.pole_width      = ${model.get(['magnet', 'magnetSectorLinear', 'pole_width'])*1e3}
m.yoke_heigth     = ${model.get(['magnet', 'magnetSectorLinear', 'yoke_height'])*1e3}
m.magn_len        = ${model.get(['magnet', 'magnetSectorLinear', 'magn_len'])*100}
m.magn_rem        = m.remanenc
m.gap_ma_yoke     = ${model.get(['magnet', 'magnetSectorLinear', 'gap_ma_yoke'])*1e3}
m.magn_ori        = ${model.get(['magnet', 'magnetSectorLinear', 'magn_ori'])}
m.airgap_shape    = ${model.get(['magnet', 'magnetSectorLinear', 'airgap_shape'])*1e3}
m.magn_type       = ${model.get(['magnet', 'magnetSectorLinear', 'magn_type'])}

m.zeroangl        =          0.000 --   Reference angle to x-axis [grad]        
 
 m.mcvkey_yoke =   mcvkey_yoke
 
 pre_models("M-Sector-Lin");
