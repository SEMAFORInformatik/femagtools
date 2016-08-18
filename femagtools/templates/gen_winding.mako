
--  Gen_winding
m.num_phases      =  ${model.get(['windings','num_phases'])}
m.num_layers      =  ${model.get(['windings','num_layers'])}
m.num_wires       =  ${model.get(['windings','num_wires'])}
m.current         =   0.0
m.coil_span       =  ${model.get(['windings','coil_span'])}

m.mat_type        =   1.0 -- rotating 
m.wind_type       =   1.0 -- winding & current
m.win_asym        =   1.0 -- sym

m.curr_inp        =   0.0 -- const
m.dq_offset       =   0

 pre_models("Gen_winding")
 pre_models("gen_pocfile") 
