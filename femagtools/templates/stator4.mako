m.yoke_diam   = dy1
m.inside_diam = da1
if( m.el_order_ag == nil ) then
  m.el_order_ag   =  1 --   El. order in air gap: lin=1: quadr=2    
end
m.slot_height     =     ${model['slot_height']*1e3}
m.slot_h1         =     ${model['slot_h1']*1e3}
m.slot_h2         =     ${model['slot_h2']*1e3}
m.slot_h3         =     ${model['slot_h3']*1e3}
m.slot_h4         =     ${model['slot_h4']*1e3}
m.slot_r1         =     ${model['slot_r1']*1e3}
m.slot_width      =     ${model['slot_width']*1e3}
m.wedge_width1    =     ${model['wedge_width1']*1e3}
m.wedge_width2    =     ${model['wedge_width2']*1e3}
m.wedge_width3    =     ${model['wedge_width3']*1e3}

m.num_layer       =     ${model['num_layers']}
m.zeroangl        =     ${model['zeroangle']}

 m.mcvkey_yoke    =    mcvkey_yoke
 m.wdg_location   =  -1.0 -- stator (internal values)

pre_models( "STATOR_4")

if mcvkey_teeth ~= nil then
  m.rlength         =     ${model.get('rlength', 1)*100}  
  if da1 > da2 then
     r = (da1 + m.slot_height)/2
  else
     r = (da1 - m.slot_height)/2
  end  
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_teeth, m.rlength)
end
