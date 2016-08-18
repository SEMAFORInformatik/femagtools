
m.magn_rad       = da2/2
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2


m.magn_height     =    ${model.get(['magnet','magnetSector', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetSector', 'magn_width_pct'])*100}
m.condshaft_r     =    ${model.get(['magnet','magnetSector', 'condshaft_r'])*1e3}
m.magn_num        =    ${model.get(['magnet','magnetSector', 'magn_num'])}
m.magn_perm       =    ${model.get(['magnet','magnetSector', 'magn_rfe'])*1e3}
m.magn_l          =    ${model.get(['magnet','magnetSector', 'magn_len'])*100}
m.magn_ori        =    ${model.get(['magnet','magnetSector', 'magn_ori'])}
m.magn_type       =    ${model.get(['magnet','magnetSector', 'magn_type'])}
m.magn_shape      =    ${model.get(['magnet','magnetSector', 'magn_shape'])*1e3}
m.br_height       =    ${model.get(['magnet','magnetSector', 'bridge_height'])*1e3}
m.br_width        =    ${model.get(['magnet','magnetSector', 'bridge_width'])*1e3}

m.zeroangl        =     0.0
m.cond_shaft      =     0.000
m.mcvkey_yoke     =     '${model.magnet.get('mcvkey_yoke','dummy')}'
m.mcvkey_mshaft   =     '${model.magnet.get('mcvkey_mshaft','dummy')}'

m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet-Sector")
