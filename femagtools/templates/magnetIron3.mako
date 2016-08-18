
m.magn_rad       = da2/2
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model.get(['magnet','magnetIron3', 'magn_height'])*1e3}
m.iron_bfe        =    ${model.get(['magnet','magnetIron3', 'iron_bfe'])*1e3}
m.gap_ma_iron     =    ${model.get(['magnet','magnetIron3', 'gap_ma_iron'])*1e3}
m.air_triangle    =    ${model.get(['magnet','magnetIron3', 'air_triangle'])}
m.iron_height     =    ${model.get(['magnet','magnetIron3', 'iron_height'])*1e3}
m.gap_ma_rigth    =    ${model.get(['magnet','magnetIron3', 'gap_ma_right'])*1e3}
m.gap_ma_left     =    ${model.get(['magnet','magnetIron3', 'gap_ma_left'])*1e3}
m.shaft_rad      =     ${model.get(['magnet','magnetIron3', 'condshaft_r'])*1e3}
m.magn_num        =    ${model.get(['magnet','magnetIron3', 'magn_num'])}
m.magn_ori        =    ${model.get(['magnet','magnetIron3', 'magn_ori'])}
m.iron_shape      =    ${model.get(['magnet','magnetIron3', 'iron_shape'])*1e3}

m.zeroangl        =     0.0

m.mcvkey_yoke     =   '${model.magnet.get('mcvkey_yoke','dummy')}'
m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet Iron 3")
