
m.l_corner_x0     = 0.0
m.l_corner_y0     = 0.0
m.slot_height     = ${model.get(['stator', 'stator3Linear', 'slot_height'])*1e3}
m.slot_h1         = ${model.get(['stator', 'stator3Linear', 'slot_h1'])*1e3}
m.slot_h2         = ${model.get(['stator', 'stator3Linear', 'slot_h2'])*1e3}
m.tip_slot        = ${model.get(['stator', 'stator3Linear', 'tip_slot'])*1e3}
m.yoke_height     = ${model.get(['stator', 'stator3Linear', 'yoke_height'])*1e3}
m.slot_r1         = ${model.get(['stator', 'stator3Linear', 'slot_r1'])*1e3}
m.slot_r2         = ${model.get(['stator', 'stator3Linear', 'slot_r2'])*1e3}
m.width_bz        = ${model.get(['stator', 'stator3Linear', 'width_bz'])*1e3}
m.tooth_width     = ${model.get(['stator', 'stator3Linear', 'tooth_width'])*1e3}

m.middle_line     = ${model.get(['stator', 'stator3Linear', 'middle_line'])}
m.zeroangl        =          0.000 --   Reference angle to x-axis [grad]        
 
 m.mcvkey_yoke = '${model.stator.get('mcvkey_yoke', 'dummy')}'                                                
 
 pre_models("STATOR3_Linear");
