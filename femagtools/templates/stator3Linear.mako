

m.cood_system     = ${model['coord_system'])} -- 1: x/y or 2: r/z 

m.l_corner_x0     = 0.0
m.l_corner_y0     = 0.0
m.slot_height     = ${model['slot_height']*1e3}
m.slot_h1         = ${model['slot_h1']*1e3}
m.slot_h2         = ${model['slot_h2']*1e3}
m.tip_slot        = ${model['tip_slot']*1e3}
m.yoke_height     = ${model['yoke_height']*1e3}
m.slot_r1         = ${model['slot_r1']*1e3}
m.slot_r2         = ${model['slot_r2']*1e3}
m.width_bz        = ${model['width_bz']*1e3}
m.tooth_width     = ${model['tooth_width']*1e3}

m.middle_line     = ${model['middle_line']}
m.zeroangl        =          0.000 --   Reference angle to x-axis [grad]        
 
 m.mcvkey_yoke = mcvkey_yoke
 
 pre_models("STATOR3_Linear");
