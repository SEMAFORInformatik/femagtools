-- calculate flux linkages, torque, and rel remanence (magstatic mode)
--
-- model:
--    Hk (kA/m) knee point field strength
--    curvec (A) phase current amplitude samples (list of [ía, ib, ic])
--    num_par_wdgs (number of parallel winding groups, default 1)
--    fc_radius (m) radius of airgap center
--
-- creates file psi-torq-rem.dat in current directory with columns:
--  displ curr1 curr2 curr3 psi1 psi2 psi3 torq rrem

<%include file="magnet-data.mako"/>

function dmg(ek, Br, alfam, murm, Bd, Bq, Hd, Hq, alfahm, Hk)
   if Hd < Hk then
      muem = murm*4*math.pi*1e-7
      Brn = Bd - Hk*1e3*muem
      if Br < 1e-5 then
         Brn = 1e-5
      end
      return Brn, alfam, 1
   end
   return Br, alfam, 1
end

function calc_flux_torq_rem(curvec)
    for k=1,3 do
      def_curr_wdg(k, curvec[k]/a, 0)
    end

    dRR = 0
    RR = 0
    maxit=m.num_nonl_it
    maxcop=m.error_perm -- err_perm in %
    permode='restore'
    repeat
       calc_field_single({
        maxit=maxit, maxcop=maxcop, -- err_perm in %
        permode=permode})
       if(maxit > 1) then
          maxit = 1
          permode='actual'
          maxcop = 0.05
       end
       stat, RR, dRR = calc_demag(5, Hk)
       --printf("%g %g", RR, dRR)
    until math.abs(dRR) < 1e-5

    psi = {}
    for k=1,3 do
      psir, psii = flux_winding_wk(k)
      psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}
    end

    fr, ft, tq, fx, fy = force_torque()

  return psi, tq, RR
end
%if model.get('fc_radius', 0):
if m.fc_radius == nil then
  m.fc_radius = ${model['fc_radius']*1e3}
end
%endif
%if type(model.get('curvec')[0]) is list:
curvec = {${','.join(['{'+','.join([str(x) for x in y])+'}' for y in model['curvec']])}} -- A
%else:
curvec = {{${','.join([str(x) for x in model['curvec']])}}} -- A
%endif
a=${model.get('num_par_wdgs', 1)}   -- parallel branches

ksym = m.num_poles/m.npols_gen

-- HcB = Brem*tempcoefbr*(magn_temp-20)+1)/muerel/12.565e-7
-- Hcmin = HcJ*tempcoefhc*(magn_temp-20.0)+1)/HcB*1e2 -- limit of demagnetization in
Hk = ${model.get('Hk', -999)}

phi = ${model.get('phi', 0)}*180/math.pi
-- initialize rotate
rotate({
    airgap = m.fc_radius,    -- air gap radius
    region = "inside",       -- region to rotate
    mode   = "save"         -- save initial model state
})

rotate({angle=phi, mode="absolute"})

file_psi = io.open("psi-torq-rem.dat","w")
for i=1, #curvec do
    psi, tq, rr = calc_flux_torq_rem(curvec[i])

    print(string.format(" current: %d/%d %g, %g, %g torque %g rr %g\n",
        i, #curvec, curvec[i][1], curvec[i][2], curvec[i][3], tq, rr))
    file_psi:write(string.format("%g ", phi))
    for k=1, 3 do
      file_psi:write(string.format("%g ", curvec[i][k]))
    end
    for k=1, 3 do
      file_psi:write(string.format("%g ", psi[k][1]))
    end
    file_psi:write(string.format("%g ", tq))
    file_psi:write(string.format("%g ", rr))
    file_psi:write("\n")

% if model.get('plots', []):
  if i == #curvec-2 then
<%include file="plots.mako"/>
  end
% endif
end

  rotate({mode = "reset"})  -- restore the initial state (discard any changes)
  file_psi:close()
