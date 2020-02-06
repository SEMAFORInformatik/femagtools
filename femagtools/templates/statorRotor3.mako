m.yoke_diam   = dy1
if (type(da1) == "table") then
  m.inside_diam = da1[2]
else
  m.inside_diam = da1
end
m.slot_height = ${model['slot_height']*1e3}
m.slot_h1     = ${model['slot_h1']*1e3}
m.slot_h2     = ${model['slot_h2']*1e3}
m.slot_width  = ${model['slot_width']*1e3}
m.slot_r1     = ${model['slot_r1']*1e3}
m.slot_r2     = ${model['slot_r2']*1e3}
m.wedge_width1= ${model['wedge_width1']*1e3}
m.wedge_width2= ${model['wedge_width2']*1e3}
m.middle_line = ${model['middle_line']}
m.tooth_width = ${model['tooth_width']*1e3}
m.slot_top_sh = ${model['slot_top_sh']}
   
m.zeroangl    = ${model['zeroangle']}
m.rlength     = ${model['rlength']*100}  

m.mcvkey_yoke = mcvkey_yoke
m.wdg_location=  -1.0 -- stator (internal values)
   
--m.delta_angle_ndchn =          0.8 --   angle nodechain airgap in [degr]
pre_models("STATOR_3")                                                      

if mcvkey_teeth ~= nil then
  if m.inside_diam > m.yoke_diam then
     r = (m.inside_diam - m.slot_height)/2
  else
     r = (m.inside_diam + m.slot_height)/2
  end  
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_teeth, m.rlength)
end
