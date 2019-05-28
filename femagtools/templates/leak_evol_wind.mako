-- End-Winding Leakage 
m.nseg            = m.tot_num_slot/m.num_slots --   Number of segments                      
m.npolsim         = m.npols_gen      --   Number of poles simulated             

m.evol1rad        = ${model['evol1rad']} --  Top radius of first evolvent [mm]  
m.evol2rad        = ${model['evol2rad']} --  Top radius of second evolvent [mm]   
m.botlevel        = ${model['botlevel']} --  Level at bottom of evolvents [mm] 
m.toplevel        = ${model['toplevel']} --  Level at top of evolvents [mm]   
m.endheight       = ${model['endheight']} --  End winding height [mm]       
m.evolbend        = ${model['evolbend']} --  Bending radius [mm] 
m.wiredia         = ${model['wiredia']} --  Wire diameter [mm]                    

 pre_models("leak_evol_wind")
