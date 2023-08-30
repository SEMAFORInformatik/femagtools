m.yoke_rad       = dy2/2
m.rotor_rad      = da2/2

m.yoke_height     =   ${model['yoke_height']*1e3}

m.iron_h1         =   ${model['iron_h1']*1e3}
m.iron_h2         =   ${model['iron_h2']*1e3}
m.iron_b          =   ${model['iron_b']*1e3}
m.magn_width      =   ${model['magn_width']*1e3}
m.magn_height     =   ${model['magn_height']*1e3}
m.iron_bfe        =   ${model['iron_bfe']*1e3}
m.iron_bfo        =   ${model['iron_bfo']*1e3}
m.iron_shape      =   ${model['iron_shape']*1e3}
m.iron_hp         =   ${model['iron_hp']*1e3}
m.magn_num        =   ${model['magn_num']}

m.nodedist        =   ${model.get('nodedist',0.0)}
m.zeroangl        =   ${model.get('zeroangle',0)}
 
m.mcvkey_yoke     =   mcvkey_yoke
 
 pre_models("Rotor_FC2")                                                     
