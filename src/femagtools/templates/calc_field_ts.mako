-- TS FE Simulation
--
m.speed = ${model['speed']*60}

-- Winding definition
----------------------
% if 'netlists' in model:
src_type = "extern"
% else:
src_type = "${model['src_type']}"
% endif
get_wdg_keys("st_wkeys")
for i=1,#st_wkeys do
  change_wdg_type(st_wkeys[i],"wire",src_type)
end

resist=0 -- ${model.get('resistance', 0)}
% if 'l_end_winding' in model:
induct=${model['l_end_winding']}
% else:
if leak ~= nil then
  induct = leak[2]
else
  induct = 0.0
end
%endif
if(m.cufilfact ~= nil) then
fillfac = 100*m.cufilfact  -- %
else
fillfac = 45
end
if(m.culength ~= nil) then
length = m.culength    -- rel. wire length in %
else
length = 140
end
if(m.cuconduct ~= nil) then
conduc = m.cuconduct
else
conduc = 56e6
end
muer = 1.0
for i = 1,#st_wkeys do
  def_ext_resist(st_wkeys[i],resist,induct)
  -- To make resultats comparable (conductivity is 0)
  def_wdg_material(st_wkeys[i], conduc, muer, fillfac, length)
end

-- Winding group definition
----------------------------
get_wdg_keys("st_wkeys")
grp_star = def_new_grp("stat",0, magenta)
for i=1,#st_wkeys do
  add_wdg_to_grp(grp_star, st_wkeys[i])
end

% if 'netlists' not in model:
-- Excitation Source Definition
-------------------------------
src_ampl = ${model['src_ampl']}
src_kind = "sinus"
% if 'freq' in model:
freq = ${model['freq']}
% else:
freq = m.speed/60*m.num_poles/2
% endif
phi = {${','.join([str(phi) for phi in model['current_angles']])}}
get_wdg_keys("st_wkeys")
for i=1,#st_wkeys do
  def_ext_source(st_wkeys[i], src_type, src_ampl, src_kind, phi[i], freq )
end
% endif

-- Start values
---------------
for i=1,#st_wkeys do
  def_curr_wdg(st_wkeys[i], 0, 0 )
end
set_speed(m.speed,0.0)

-- Calculation control
----------------------
inertia = 0.01
m0 = 0.0
m1 = 0.0
m2 = 0.0
m3 = 0.0
m4 = 0.0
set_ext_torque(inertia, m0, m1,m2, m3, m4)

sl_interpolation = 0
speed_control = 0
rho_0 = 0.0
beta_0 = 0.0
set_calc_mode(speed_control, rho_0, beta_0)

-- FE Simulation
----------------
cond_temp = ${model.get('cond_temp', 20)}
set_dev_data("cond_temp", cond_temp, cond_temp)
m.magn_temp = ${model.get('magn_temp', 20)}
set_dev_data("magn_temp", m.magn_temp)

% if 'netlists' in model:
function write_netlist(branches)
   nlf = io.open(model .. '.net', 'w')
   nlf:write('{"Branches":[\n')
   for i=1, #branches-1, 1 do
      nlf:write(string.format("%s,\n", branches[i]))
   end
   if #branches > 0
   then
      nlf:write(string.format("%s\n", branches[#branches]))
   end
   nlf:write("]}")
   io.close(nlf)
end

<% import json %>
branches = {
${',\n'.join([f"\'{json.dumps(b)}\'" for b in model['netlists'][0]['Branches']])}
}
write_netlist(branches)
% endif
% if model.get('vtu_movie', False):
set_store_mode("stx+vtu")
% else:
set_store_mode("stx")
% endif
% if isinstance(model['sim_time'], list):
sim_time = ${model['sim_time'][0]}
% else:
sim_time = ${model['sim_time']}
% endif
store_time = ${model.get('store_time', 0)}
dtmin = ${model['dtmin']}
dtmax = ${model['dtmax']}
resmin = ${model.get('resmin', 1e-5)}
resmax = ${model.get('resmax', 1e-4)}
mode = 0
t_end = calc_field_ts(mode, sim_time, store_time, dtmin, dtmax, resmin, resmax)
% if isinstance(model['sim_time'], list):
sim_time = ${model['sim_time'][1]}

% if 'netlists' in model and len(model['netlists'])>1:
branches = {
${',\n'.join([f"\'{json.dumps(b)}\'" for b in model['netlists'][1]['Branches']])}
}
write_netlist(branches)
% endif
mode = 4
t_end = calc_field_ts(mode, sim_time, store_time, dtmin, dtmax, resmin, resmax)
% endif
