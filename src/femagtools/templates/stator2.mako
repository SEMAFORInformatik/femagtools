m.yoke_rad   = dy1/2
m.inside_rad = da1/2

m.slot_t1         =     ${model['slot_t1']*1e3}
m.slot_t2         =     ${model['slot_t2']*1e3}
m.slot_t3         =     ${model['slot_t3']*1e3}
m.slot_depth      =     ${model['slot_depth']*1e3}
m.slot_width      =     ${model['slot_width']*1e3}
m.corner_width    =     ${model['corner_width']*1e3}
m.num_layer       =     ${model['num_layers']}
m.zeroangl        =     ${model['zeroangle']}
m.rlength        =      ${model['rlength']*100}

 m.mcvkey_yoke    =    mcvkey_yoke
 m.wdg_location   =  -1.0 -- stator (internal values)

pre_models( "STATOR_2")

