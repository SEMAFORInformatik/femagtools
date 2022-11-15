-- End-Winding Leakage
m.nseg            = m.tot_num_slot/m.num_slots --   Number of segments
m.npolsim         = m.npols_gen      --   Number of poles simulated
% if 'perimrad' not in model:
m.perimrad        = (da1+dy1)/4   --   Radius of perimeter [mm]
% else:
m.perimrad        = ${model['perimrad']*1e3}   --   Radius of perimeter [mm]
% endif
m.vbendrad        = ${model['vbendrad']*1e3} --   Bending radius vertical [mm]
m.endheight       = ${model['endheight']*1e3} --   End winding height [mm]
m.wiredia         = ${model['wiredia']*1e3} --   Wire diameter [mm]
pre_models("leak_dist_wind")

 post_models("end_wind_leak","leak")
file_leak = io.open("end_wind_leak.dat","w")
file_leak:write(string.format("%g %g %g\n",
                leak[1], leak[2], leak[3]))
file_leak:close()
