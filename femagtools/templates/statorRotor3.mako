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

m.mcvkey_yoke = '${model.stator.get('mcvkey_yoke','dummy')}'  
m.wdg_location=  -1.0 -- stator (internal values)
   
--m.delta_angle_ndchn =          0.8 --   angle nodechain airgap in [degr]
pre_models("STATOR_3")                                                      
