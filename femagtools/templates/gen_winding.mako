
--  Gen_winding
m.num_phases      =  ${model.get(['windings','num_phases'])}
m.num_layers      =  ${model.get(['windings','num_layers'])}
m.num_wires       =  ${model.get(['windings','num_wires'])}
m.coil_span       =  ${model.get(['windings','coil_span'])}

% if model.get('move_action') == 0:
m.current         =   0.0
m.mat_type        =   1.0 -- rotating 
m.wind_type       =   1.0 -- winding & current
m.win_asym        =   1.0 -- sym

m.curr_inp        =   0.0 -- const
m.dq_offset       =   0

 pre_models("Gen_winding")
 pre_models("gen_pocfile") 
% else:
function get_coil(tooth, slots, pols)
  angle = ((tooth-1)*180.0*pols/slots)%360

  if (angle >= 330 or angle <  30) then 
    return 1
  end
  if (angle >= 30 and angle < 90) then
    return -3 
  end
  if (angle >= 90 and angle < 150) then 
    return 2
  end
  if (angle >= 150 and angle < 210) then
    return -1 
  end
  if (angle >= 210 and angle < 270) then 
    return 3
  end
  if (angle >= 270 and angle < 330) then
    return -2
  end
end

xl=-m.width_bz/3.0
xr=m.width_bz/3.0
ys =10.0
color={green, yellow, magenta}
wkey={0,0,0}

for i=1,m.num_phases do
  for j=1,m.tot_num_slot do
    xsl = xl + (j-1)*m.width_bz
    xsr = xr + (j-1)*m.width_bz
    coil = get_coil(j, m.tot_num_slot, m.num_poles)
    if (i==math.abs(coil)) then
      wdgid = "w"..i
      if (coil > 0) then
        dirl = wi
        dirr = wo
      else
        dirl = wo
        dirr = wi
      end
  
      if (xsl > 0 and xsl < m.num_sl_gen*m.width_bz) then
        if (wkey[i] == 0) then
          wkey[i]=def_new_wdg(xsl,ys, color[i], wdgid, m.num_wires, 0, dirl)
        else
          add_to_wdg(xsl,ys, wkey[i], dirl, wser)
        end
      end
      if (xsr > 0 and xsr < m.num_sl_gen*m.width_bz) then
        if (wkey[i] == 0) then
          wkey[i]=def_new_wdg(xsr,ys, color[i], wdgid, m.num_wires, 0, dirr)
        else
          add_to_wdg(xsr,ys, wkey[i], dirr, wser)
        end
      end
    end
  end
end


----------------------------------------------------------------------
-- create poc-File 
----------------------------------------------------------------------
phases = 3
base_angle = 0.0
period = 2*m.tot_num_slot*m.width_bz/m.num_poles
f = assert(io.open(model..'_'..m.num_poles.."p.poc","w"));
  f:write(string.format("          %i\n",phases));
  for i = 1, phases do
    if i<10      then f:write(string.format("          %i\n",i));
    elseif i<100 then f:write(string.format("         %i\n",i));
    end
  end
  for i = 0,phases-1 do
    angle = base_angle+i*360.0/phases
    if     angle<10   then f:write(string.format("   %10.8f\n",angle));
    elseif angle<100  then f:write(string.format("   %10.7f\n",angle));
    elseif angle<1000 then f:write(string.format("   %10.6f\n",angle));
    end
  end
  if     period<10   then f:write(string.format("   %10.8f\n",period));
  elseif period<100  then f:write(string.format("   %10.7f\n",period));
  elseif period<1000 then f:write(string.format("   %10.6f\n",period));
  end
  f:write("sin\n");
  f:write("   0.00000000\n");
  f:write("          0\n");
io.close(f);

%endif
