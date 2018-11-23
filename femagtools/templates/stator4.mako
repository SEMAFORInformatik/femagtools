m.yoke_diam   = dy1
m.inside_diam = da1

m.slot_height     =     ${model.get(['stator', 'stator4', 'slot_height'])*1e3}
m.slot_h1         =     ${model.get(['stator', 'stator4', 'slot_h1'])*1e3}
m.slot_h2         =     ${model.get(['stator', 'stator4', 'slot_h2'])*1e3}
m.slot_h3         =     ${model.get(['stator', 'stator4', 'slot_h3'])*1e3}
m.slot_h4         =     ${model.get(['stator', 'stator4', 'slot_h4'])*1e3}
m.slot_r1         =     ${model.get(['stator', 'stator4', 'slot_r1'])*1e3}
m.slot_width      =     ${model.get(['stator', 'stator4', 'slot_width'])*1e3}
m.wedge_width1    =     ${model.get(['stator', 'stator4', 'wedge_width1'])*1e3}
m.wedge_width2    =     ${model.get(['stator', 'stator4', 'wedge_width2'])*1e3}
m.wedge_width3    =     ${model.get(['stator', 'stator4', 'wedge_width3'])*1e3}

m.num_layer       =     ${model.get(['windings', 'num_layers'])}
m.zeroangl        =     ${model.stator.get('zeroangle', 0.0)}

 m.mcvkey_yoke    =    mcvkey_yoke
 m.wdg_location   =  -1.0 -- stator (internal values)

pre_models( "STATOR_4")

if mcvkey_teeth ~= nil then
  if da1 > da2 then
     r = (da1 + m.slot_height)/2
  else
     r = (da1 - m.slot_height)/2
  end  
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, blue, mcvkey_teeth, m.rlength)
end
