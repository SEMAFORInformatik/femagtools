-- calculate noload flux and airgap flux density
--
-- model:
--    currvec, voltvec (A rms, V rms
--    num_par_wdgs (number of parallel winding groups, default 1)
--
function calc_field(curvec, psifile, bagfile, phiagfile)
  psi={}
  cur={}
  volt={}
  bag={}
  phiag = {}
  f1=0
  a=${model.get('num_par_wdgs', 1)}   -- parallel branches

  file_psi = io.open(psifile,"w")

  for i = 1, #curvec do
    amp = math.sqrt(2)*curvec[i]/a
    for k=1,3 do
    --Set currents
      phi = (k-1)*2*math.pi/3
      cur[k] = {amp*math.cos(phi), 0}
--                 -amp*math.sin(phi)}
--      def_volt_wdg(k,volt[k][1],volt[k][2])
      def_curr_wdg(k, cur[k][1], cur[k][2])
    end
    calc_field_single({maxit=maxit, maxcop=du_u0,
        permode='restore',freq=f1})

    for k=1, 3 do
      psir,psii=flux_winding_wk(k)
      curr, curi = get_wdg_data("cur", k)
      voltr, volti = get_wdg_data("volt", k)
      cur[k] = {a*curr, a*curi}
      volt[k] = {ksym*m.arm_length*voltr/a, ksym*m.arm_length*volti/a}
      psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}

--[[
      math.sqrt(volt[1][1]^2 + volt[1][2]^2),
      math.sqrt(volt[2][1]^2 + volt[2][2]^2),
      math.sqrt(volt[3][1]^2 + volt[3][2]^2),
--]]

    end
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
    --effval,
        math.sqrt((cur[1][1]^2 + cur[1][2]^2)/2),
        math.sqrt((cur[2][1]^2 + cur[2][2]^2)/2),
        math.sqrt((cur[3][1]^2 + cur[3][2]^2)/2),
        math.sqrt((psi[1][1]^2 + psi[1][2]^2)/2),
        math.sqrt((psi[2][1]^2 + psi[2][2]^2)/2),
        math.sqrt((psi[3][1]^2 + psi[3][2]^2)/2)))
--[[
    m.coord_x1, m.coord_y1 = pd2c( m.fc_radius, 0)
    m.coord_x2, m.coord_y2 = pd2c( m.fc_radius, 360/m.num_poles)
    post_models("flux(x1,x2)","bflux")
      phiag[i] = {ksym*bflux[2]/a*m.arm_length, ksym*bflux[2]/a*m.arm_length}
--]]
    post_models("induct(x)","b")    -- Calculate field distribution
      bag[i] = {}
      n = 1
      for j = 1, table.getn(b), 3 do
       bag[i][n] = {b[j],b[j+1],b[j+2]}
       n = n+1
      end
  end
  file_psi:close()

  file_bag = io.open(bagfile,"w")
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
--[[
  file_phiag = io.open(phiagfile,"w")
  for i=1, #phiag do
      file_phiag:write(string.format("%g %g\n",
        phiag[i][1], phiag[i][2]))
  end
  file_phiag:close()
  --]]
end

ksym = m.num_poles/m.npols_gen
maxit = m.num_nonl_it
du_u0 = m.error_perm
rmod = m.perm_mode
cmod = 0
print("\n No load flux simulation (DC)\n")

curvec = {${','.join([str(x) for x in model['curvec']])}} -- A rms

calc_field(curvec, "noloadflux.dat", "noloadbag.dat", "psiairgap.flux")
