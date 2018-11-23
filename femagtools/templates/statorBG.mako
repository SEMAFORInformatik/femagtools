 
m.yoke_diam       = dy1
m.inside_diam     = da1
m.yoke_diam_ins   = ${model.get(['stator','statorBG','yoke_diam_ins'])*1e3}
m.slot_h1         = ${model.get(['stator','statorBG','slot_h1'])*1e3}
m.slot_h3         = ${model.get(['stator','statorBG','slot_h3'])*1e3}
m.slot_width      = ${model.get(['stator','statorBG','slot_width'])*1e3}
m.slot_r1         = ${model.get(['stator','statorBG','slot_r1'])*1e3}
m.slot_r2         = ${model.get(['stator','statorBG','slot_r2'])*1e3}
m.tooth_width     = ${model.get(['stator','statorBG','tooth_width'])*1e3}
m.middle_line     = ${model.get(['stator','statorBG','middle_line'])}

m.tip_rad         = ${model.get(['stator','statorBG','tip_rad'])*1e3}
m.slottooth       = ${model.get(['stator','statorBG','slottooth'])*1e3}
m.wdg_location    =  -1.0 -- stator (internal values)
 
m.zeroangl    = ${model.stator.get('zeroangle',0)}

m.rlength     = ${model.stator.get('rlength',1.0)*100}  
m.mcvkey_yoke = mcvkey_yoke
 
 pre_models("STATOR_BG")

if mcvkey_teeth ~= nil then
  r = (da1 + m.yoke_diam_ins)/4
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, blue, mcvkey_teeth, m.rlength)
end
