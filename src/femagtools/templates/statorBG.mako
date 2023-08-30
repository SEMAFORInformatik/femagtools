 
m.yoke_diam       = dy1
m.inside_diam     = da1
m.yoke_diam_ins   = ${model['yoke_diam_ins']*1e3}
m.slot_h1         = ${model['slot_h1']*1e3}
m.slot_h3         = ${model['slot_h3']*1e3}
m.slot_width      = ${model['slot_width']*1e3}
m.slot_r1         = ${model['slot_r1']*1e3}
m.slot_r2         = ${model['slot_r2']*1e3}
m.tooth_width     = ${model['tooth_width']*1e3}
m.middle_line     = ${model['middle_line']}

m.tip_rad         = ${model['tip_rad']*1e3}
m.slottooth       = ${model['slottooth']*1e3}
m.wdg_location    =  -1.0 -- stator (internal values)
 
m.zeroangl    = ${model['zeroangle']}

m.rlength     = ${model['rlength']*100}  
m.mcvkey_yoke = mcvkey_yoke
 
 pre_models("STATOR_BG")

if mcvkey_teeth ~= nil then
  r = (da1 + m.yoke_diam_ins)/4
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_teeth, m.rlength)
end
