m.yoke_diam   = dy1
m.inside_diam = da1
m.slot_height = ${model.get(['stator','statorRotor3','slot_height'])*1e3}
m.slot_h1     = ${model.get(['stator','statorRotor3','slot_h1'])*1e3}
m.slot_h2     = ${model.get(['stator','statorRotor3','slot_h2'])*1e3}
m.slot_width  = ${model.get(['stator','statorRotor3','slot_width'])*1e3}
m.slot_r1     = ${model.get(['stator','statorRotor3','slot_r1'])*1e3}
m.slot_r2     = ${model.get(['stator','statorRotor3','slot_r2'])*1e3}
m.wedge_width1= ${model.get(['stator','statorRotor3','wedge_width1'])*1e3}
m.wedge_width2= ${model.get(['stator','statorRotor3','wedge_width2'])*1e3}
m.middle_line = ${model.get(['stator','statorRotor3','middle_line'])}
m.tooth_width = ${model.get(['stator','statorRotor3','tooth_width'])*1e3}
m.slot_top_sh = ${model.get(['stator','statorRotor3','slot_top_sh'])}
   
m.zeroangl    = ${model.stator.get('zeroangle',0)}
m.rlength     = ${model.stator.get('rlength',1.0)*100}  

m.mcvkey_yoke = mcvkey_yoke
m.wdg_location=  -1.0 -- stator (internal values)
   
--m.delta_angle_ndchn =          0.8 --   angle nodechain airgap in [degr]
pre_models("STATOR_3")                                                      

if mcvkey_teeth ~= nil then
  if da1 > da2 then
     r = (da1 + m.slot_height)/2
  else
     r = (da1 - m.slot_height)/2
  end  
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, blue, mcvkey_teeth, m.rlength)
end
