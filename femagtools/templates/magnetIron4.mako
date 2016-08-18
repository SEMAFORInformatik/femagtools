
m.magn_rad       = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model.get(['magnet','magnetIron4', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetIron4', 'magn_width'])*1e3}
m.gap_ma_iron     =    ${model.get(['magnet','magnetIron4', 'gap_ma_iron'])*1e3}
m.iron_shape      =    ${model.get(['magnet','magnetIron4', 'iron_shape'])*1e3}
m.air_space_h     =    ${model.get(['magnet','magnetIron4', 'air_space_h'])*1e3}
m.iron_bfe        =    ${model.get(['magnet','magnetIron4', 'iron_bfe'])*1e3}
m.magn_di_ra      =    ${model.get(['magnet','magnetIron4', 'magn_di_ra'])*1e3}
m.corner_r        =    ${model.get(['magnet','magnetIron4', 'corner_r'])*1e3}
m.air_sp_ori      =    ${model.get(['magnet','magnetIron4', 'air_sp_ori'])}
m.magn_ori        =    ${model.get(['magnet','magnetIron4', 'magn_ori'])}
m.magn_num        =    ${model.get(['magnet','magnetIron4', 'magn_num'])}
m.cond_shaft      =    0.0 -- ignored
m.zeroangl        =     0.0

m.mcvkey_yoke     =   '${model.magnet.get('mcvkey_yoke','dummy')}'
m.mcvkey_mshaft    =   '${model.magnet.get('mcvkey_mshaft','dummy')}'
m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet Iron 4")
