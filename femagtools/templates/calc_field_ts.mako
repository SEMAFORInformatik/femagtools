-- TS FE Simulation
--
-- Winding definition
----------------------
--Statorwicklungs-Keys als extern definieren
src_type = "${model['src_type']}"
get_wdg_keys("st_wkeys")
for i=1,#st_wkeys do
  change_wdg_type(st_wkeys[i],"wire",src_type)
end

resist=0.0
induct=0.0
for i = 1,#st_wkeys do
  def_ext_resist(st_wkeys[i],resist,induct)
  -- Damit Resultate vergleichbar sind
  -- FEMAG-DC setzt bei Berechnung die Leitfähigkeit der Wicklung auf 0
  conduc=56e6
  --conduc=0
  muer=1.0
  fillfac=45
  length=180
  def_wdg_material(st_wkeys[i], conduc, muer, fillfac, length)
end

-- Winding group definition
----------------------------
get_wdg_keys("st_wkeys")
grp_star = def_new_grp("stat",0, magenta)
for i=1,#st_wkeys do
  add_wdg_to_grp(grp_star, st_wkeys[i])
end

-- Excitation Source Definition
-------------------------------
src_ampl = ${model['src_ampl']}
src_kind = "sinus"
freq = ${model['freq']}
phi = {60,180,300}
get_wdg_keys("st_wkeys")
for i=1,#st_wkeys do
  def_ext_source(st_wkeys[i], src_type, src_ampl, src_kind, phi[i], freq )
end

-- Start values
---------------
for i=1,#st_wkeys do
  def_curr_wdg(st_wkeys[i], 0, 0 )
end
set_speed(${model['speed']*60},0.0)

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

set_store_mode("stx")

-- FE Simulation
----------------
sim_time = {${model['sim_time'][0]}, ${model['sim_time'][1]}}

store_time = ${model['store_time']}
dtmin = ${model['dtmin']}
dtmax = ${model['dtmax']}
resmin = ${model['resmin']}
resmax = ${model['resmax']}
mode = 0
t_end = calc_field_ts(mode, sim_time[1], store_time, dtmin, dtmax, resmin, resmax)

set_store_mode("vtu")
mode = 3
t_end = calc_field_ts(mode, sim_time[2], store_time, dtmin, dtmax, resmin, resmax)
