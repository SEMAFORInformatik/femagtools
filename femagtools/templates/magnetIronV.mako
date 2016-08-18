
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model.get(['magnet','magnetIronV', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetIronV', 'magn_width'])*1e3}
m.magn_angle      =    ${model.get(['magnet','magnetIronV', 'magn_angle'])}
m.magn_num        =    ${model.get(['magnet','magnetIronV', 'magn_num'])}
m.iron_hs         =    ${model.get(['magnet','magnetIronV', 'iron_hs'])*1e3}
m.iron_height     =    ${model.get(['magnet','magnetIronV', 'iron_height'])*1e3}
m.iron_shape      =    ${model.get(['magnet','magnetIronV', 'iron_shape'])*1e3}
m.air_triangle    =    ${model.get(['magnet','magnetIronV', 'air_triangle'])}
m.gap_ma_iron     =    ${model.get(['magnet','magnetIronV', 'gap_ma_iron'])*1e3}
m.magn_rem         =    ${model.get(['magnet','magnetIronV', 'magn_rem'])}
m.shaft_rad       =     ${model.get(['magnet','magnetIronV', 'condshaft_r'])*1e3}

m.zeroangl        =     0.0

m.mcvkey_yoke     =     '${model.magnet.get('mcvkey_yoke','dummy')}'

m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet Iron V")





