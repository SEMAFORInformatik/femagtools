m.yoke_rad       = dy2/2
m.rotor_rad      = da2/2

m.yoke_height     =   ${model.get(['magnet','magnetFC2', 'yoke_height'])*1e3}

m.iron_h1         =   ${model.get(['magnet','magnetFC2', 'iron_h1'])*1e3}
m.iron_h2         =   ${model.get(['magnet','magnetFC2', 'iron_h2'])*1e3}
m.iron_b          =   ${model.get(['magnet','magnetFC2', 'iron_b'])*1e3}
m.magn_width      =   ${model.get(['magnet','magnetFC2', 'magn_width'])*1e3}
m.magn_height     =   ${model.get(['magnet','magnetFC2', 'magn_height'])*1e3}
m.iron_bfe        =   ${model.get(['magnet','magnetFC2', 'iron_bfe'])*1e3}
m.iron_bfo        =   ${model.get(['magnet','magnetFC2', 'iron_bfo'])*1e3}
m.iron_shape      =   ${model.get(['magnet','magnetFC2', 'iron_shape'])*1e3}
m.iron_hp         =   ${model.get(['magnet','magnetFC2', 'iron_hp'])*1e3}
m.magn_num        =   ${model.get(['magnet','magnetFC2', 'magn_num'])}

m.nodedist        =   ${model.magnet.get('nodedist',0.0)}
m.zeroangl        =          0.000 --   Reference angle to x-axis [grad]        
 
m.mcvkey_yoke     =   '${model.magnet.get('mcvkey_yoke','dummy')}'
 
 pre_models("Rotor_FC2")                                                     
