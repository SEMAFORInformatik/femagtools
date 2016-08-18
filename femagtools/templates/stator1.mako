m.yoke_rad   = dy1/2
m.inside_rad = da1/2

m.tip_rh2         =     ${model.get(['stator', 'stator1', 'tip_rh2'])*1e3}
m.tip_rh1         =     ${model.get(['stator', 'stator1', 'tip_rh1'])*1e3}
m.slot_rf1        =     ${model.get(['stator', 'stator1', 'slot_rf1'])*1e3}
m.slot_rf2        =     m.slot_rf1
m.yoke_rad2       =     m.yoke_rad
m.slot_width      =     ${model.get(['stator', 'stator1', 'slot_width'])*1e3}
m.tooth_width     =     ${model.get(['stator', 'stator1', 'tooth_width'])*1e3}
m.zeroangl        =     ${model.stator.get('zeroangle',0.0)}

 m.mcvkey_yoke =   '${model.stator.get('mcvkey_yoke', 'dummy')}'  
 m.wdg_location    =  -1.0 -- stator (internal values)

pre_models( "STATOR_1")
