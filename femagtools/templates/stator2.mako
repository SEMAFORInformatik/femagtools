m.yoke_rad   = dy1/2
m.inside_rad = da1/2

m.slot_t1         =     ${model.get(['stator', 'stator2', 'slot_t1'])*1e3}
m.slot_t2         =     ${model.get(['stator', 'stator2', 'slot_t2'])*1e3}
m.slot_t3         =     ${model.get(['stator', 'stator2', 'slot_t3'])*1e3}
m.slot_depth      =     ${model.get(['stator', 'stator2', 'slot_depth'])*1e3}
m.slot_width      =     ${model.get(['stator', 'stator2', 'slot_width'])*1e3}
m.corner_width    =     ${model.get(['stator', 'stator2', 'corner_width'])*1e3}
m.num_layer       =     ${model.get(['windings', 'num_layers'])}
m.zeroangl        =     ${model.stator.get('zeroangle',0.0)}

 m.mcvkey_yoke    =    mcvkey_yoke
 m.wdg_location   =  -1.0 -- stator (internal values)

pre_models( "STATOR_2")

