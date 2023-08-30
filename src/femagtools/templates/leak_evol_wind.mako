-- End-Winding Leakage 
m.nseg            = m.tot_num_slot/m.num_slots --   Number of segments                      
m.npolsim         = m.npols_gen      --   Number of poles simulated             

m.evol1rad        = ${model['evol1rad']*1e3} --  Top radius of first evolvent [mm]
m.evol2rad        = ${model['evol2rad']*1e3} --  Top radius of second evolvent [mm]
m.botlevel        = ${model['botlevel']*1e3} --  Level at bottom of evolvents [mm]
m.toplevel        = ${model['toplevel']*1e3} --  Level at top of evolvents [mm]
m.endheight       = ${model['endheight']*1e3} --  End winding height [mm]
m.evolbend        = ${model['evolbend']*1e3} --  Bending radius [mm]
m.wiredia         = ${model['wiredia']*1e3} --  Wire diameter [mm]

 pre_models("leak_evol_wind")
