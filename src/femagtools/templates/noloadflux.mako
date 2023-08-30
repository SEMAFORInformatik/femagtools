-- calculate noload flux and airgap flux density
--
-- model:
--    frequency (Hz),
--    currvec, voltvec (A rms, V rms
--    num_par_wdgs (number of parallel winding groups, default 1)
--
ksym = m.num_poles/m.npols_gen
maxit = m.num_nonl_it
du_u0 = m.error_perm
rmod = m.perm_mode
cmod = 0
print("\n No load flux simulation\n")
f1 = ${model.get('frequency', 0)} -- frequency
% if model.get('frequency', 0):
  state_of_problem('mag_dynamic')
  -- set conductivity of rotor bars to 0 for noload flux
  get_sreg_keys("srkeys")
  for j=1, #srkeys do
    srname = get_sreg_data('name', srkeys[j])
    if srname == 'Bar ' then
      sekeys=get_sreg_data("sekeys",srkeys[j])
      ndkeys=get_spel_data("ndkeys",sekeys[1])
      xm, ym = 0, 0
      for i=1, #ndkeys do
        x, y = get_node_data("xy", ndkeys[i])
        xm = xm + x
        ym = ym + y
      end
      x,y=xm/#ndkeys, ym/#ndkeys
      set_mat_cond(x,y, 0, 100.0 )
    end
  end
% endif
psi={}
cur={}
volt={}
bag = {}
a=${model.get('num_par_wdgs', 1)}   -- parallel branches

file_psi = io.open("noloadflux.dat","w")
file_ind = io.open("inductances.dat","w")

% if model.get('voltvec', []):
get_wdg_keys("wdgkeys") -- get all winding keys
for j=1, #wdgkeys do
  change_wdg_type(wdgkeys[j], "wire", "voltage")
end
voltvec = {${','.join([str(x) for x in model['voltvec']])}}
for i = 1, #voltvec do
  effval = voltvec[i]
  amp = a*math.sqrt(2)*voltvec[i]/ksym/m.arm_length*1e3
  for k=1,3 do
  --Set voltages and clear currents
      phi = (k-1)*2*math.pi/3
      volt[k] = {amp*math.cos(phi),
                 -amp*math.sin(phi)}
      def_volt_wdg(k,volt[k][1],volt[k][2])
      def_curr_wdg(k, 0, 0)
  end
  calc_field_single({maxit=maxit, maxcop=du_u0,
        permode='restore',freq=f1})

  for k=1, 3 do
    psir,psii=flux_winding_wk(k)
    curr, curi = get_wdg_data("cur", k)
    cur[k] = {a*curr, a*curi}
    psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}

--[[
      math.sqrt(volt[1][1]^2 + volt[1][2]^2),
      math.sqrt(volt[2][1]^2 + volt[2][2]^2),
      math.sqrt(volt[3][1]^2 + volt[3][2]^2),
--]]

  end
% else:
curvec = {${','.join([str(x) for x in model['curvec']])}}
for i = 1, #curvec do
  amp = math.sqrt(2)*curvec[i]/a
  for k=1,3 do
  --Set currents
      phi = (k-1)*2*math.pi/3
      ire, iim = amp*math.cos(phi), amp*math.sin(phi)
      if(f1 > 0) then
         def_curr_wdg(k, ire, iim)
      else
         def_curr_wdg(k, ire)
      end
  end
  calc_field_single({maxit=maxit, maxcop=du_u0,
        permode='restore',freq=f1})

  for k=1, 3 do
    psir,psii=flux_winding_wk(k)
    voltr, volti = get_wdg_data("volt", k)
    volt[k] = {voltr, volti}
    curr, curi = get_wdg_data("cur", k)
    cur[k] = {a*curr, a*curi}
    psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}

--[[
      math.sqrt(volt[1][1]^2 + volt[1][2]^2),
      math.sqrt(volt[2][1]^2 + volt[2][2]^2),
      math.sqrt(volt[3][1]^2 + volt[3][2]^2),
--]]
  end
% endif
  for k=1, 3 do
    file_psi:write(string.format("%g %g ",
      cur[k][1], cur[k][2]))
  end
  for k=1, 3 do
    file_psi:write(string.format("%g %g ",
      volt[k][1], volt[k][2]))
  end
  for k=1, 3 do
    file_psi:write(string.format("%g %g ",
      psi[k][1], psi[k][2]))
  end
  file_psi:write("\n")
    print(string.format("curr: %g, %g, %g flux: %g, %g, %g",
      math.sqrt((cur[1][1]^2 + cur[1][2]^2)/2),
      math.sqrt((cur[2][1]^2 + cur[2][2]^2)/2),
      math.sqrt((cur[3][1]^2 + cur[3][2]^2)/2),
      math.sqrt((psi[1][1]^2 + psi[1][2]^2)/2),
      math.sqrt((psi[2][1]^2 + psi[2][2]^2)/2),
      math.sqrt((psi[3][1]^2 + psi[3][2]^2)/2)))

    post_models("induct(x)","b")    -- Calculate field distribution
      bag[i] = {}
      n = 1
      for j = 1, table.getn(b), 3 do
       bag[i][n] = {b[j],b[j+1],b[j+2]}
       n = n+1
      end

  for k=1,3 do -- self and mutual inductances
  --Set currents
     if(k==1) then
      ire, iim = 1, 0
     else
      ire, iim = 0, 0
     end
      def_curr_wdg(k, ire, iim)
  end
  calc_field_single({maxit=1, permode='actual'})

  for k=1, 3 do
    psir,psii=flux_winding_wk(k)
    psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}
  end
    -- write inductances
    file_ind:write(string.format("%g %g %g\n",
      math.sqrt((psi[1][1]^2 + psi[1][2]^2)),
      math.sqrt((psi[2][1]^2 + psi[2][2]^2)),
      math.sqrt((psi[3][1]^2 + psi[3][2]^2))))

end
file_psi:close()
file_ind:close()

file_bag = io.open("noloadbag.dat","w")
for i=1, n-1 do
  file_bag:write(string.format("%g ",
       bag[1][i][1]))
  for j=1, #bag do
    file_bag:write(string.format("%g ",
       bag[j][i][2]))
  end
  file_bag:write("\n")
end
file_bag:close()
