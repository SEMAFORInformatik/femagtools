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

mu0 = 4*math.pi*1e-7

<%include file="demagmod/dmg.lua"/>

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
        if( Hk < 0) then
          if(maxit > 1) then
            maxit = 1
            permode='actual'
            maxcop = 0.05
         end
         stat, RR, dRR = calc_demag(5, Hk)
         --printf("%g %g", RR, dRR)
       end
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
%if 'phi' in model:
phi = ${model.get('phi', 0)}*180/math.pi
-- initialize rotate
rotate({
    airgap = m.fc_radius,    -- air gap radius
    region = "inside",       -- region to rotate
    mode   = "save"         -- save initial model state
})
rotate({angle=phi, mode="absolute"})
%endif
%else:
curvec = {${','.join([str(x) for x in model['curvec']])}} -- A
phi = 90
curr = {}
curr_angles = {${','.join([str(phi) for phi in model['current_angles']])}}
%endif

a=${model.get('num_par_wdgs', 1)}   -- parallel branches

ksym = m.num_poles/m.npols_gen

-- HcB = Brem*tempcoefbr*(magn_temp-20)+1)/muerel/12.565e-7
-- Hcmin = HcJ*tempcoefhc*(magn_temp-20.0)+1)/HcB*1e2 -- limit of demagnetization in
Hk = ${model.get('Hk', -999)}

file_psi = io.open("psi-torq-rem.dat","w")
for i=1, #curvec do
  if( curr_angles ~= nil ) then
     for k=1, 3 do
        curr[k] = curvec[i]*math.sin(math.pi*(
                                         phi-curr_angles[k])/180)
     end
  else
     curr = curvec[i]
  end
  psi, tq, rr = calc_flux_torq_rem(curr)

    print(string.format(" current: %d/%d %g, %g, %g torque %g rr %g",
        i, #curvec, curr[1], curr[2], curr[3], tq, rr))

    file_psi:write(string.format("%g ", phi))
    for k=1, 3 do
      file_psi:write(string.format("%g ", curr[k]))
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
%if 'phi' in model:
rotate({mode = "reset"})
%endif
file_psi:close()
