-- End-Winding Leakage 
m.nseg            = m.tot_num_slot/m.num_slots --   Number of segments                      
m.npolsim         = m.npols_gen      --   Number of poles simulated

m.endheight       = ${model['endheight']*1e3}  -- End winding height [mm]  
m.bendrad         = ${model['bendrad']*1e3}  -- Bending radius [mm]  
m.wiredia         = ${model['wiredia']*1e3} -- Wire diameter [mm] 

 pre_models("leak_tooth_wind")
