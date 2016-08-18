
m.magn_rad       = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model.get(['magnet','magnetIron', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetIron', 'magn_width'])*1e3}
m.gap_ma_iron     =    ${model.get(['magnet','magnetIron', 'gap_ma_iron'])*1e3}
m.air_triangle    =    ${model.get(['magnet','magnetIron', 'air_triangle'])}
m.iron_height     =    ${model.get(['magnet','magnetIron', 'iron_height'])*1e3}
m.magn_rem        =    ${model.get(['magnet','magnetIron', 'magn_rem'])}
m.shaft_rad      =     ${model.get(['magnet','magnetIron', 'condshaft_r'])*1e3}
m.magn_ori        =    ${model.get(['magnet','magnetIron', 'magn_ori'])}
m.bridge_height   =    ${model.get(['magnet','magnetIron', 'bridge_height'])*1e3}
m.bridge_width    =    ${model.get(['magnet','magnetIron', 'bridge_width'])*1e3}
m.iron_shape      =    ${model.get(['magnet','magnetIron', 'iron_shape'])*1e3}
m.magn_ori        =    ${model.get(['magnet','magnetIron', 'magn_ori'])}

m.zeroangl        =     0.0

m.mcvkey_yoke     =   '${model.magnet.get('mcvkey_yoke','dummy')}'
m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet in Iron")
