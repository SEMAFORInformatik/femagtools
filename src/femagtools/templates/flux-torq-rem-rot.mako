-- calculate flux, torque, remance (magstatic mode)
--
-- model:
--    currvec (A rms)
--    num_par_wdgs (number of parallel winding groups, default 1)
--
function gcd(a, b)
	return b==0 and a or gcd(b,a%b)
end

%if model.get('sim_demagn', 0):
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
%endif

function calc_flux_torq(phi, curvec)
  psivec={}
  tqvec={}
  rrvec={}
  for i=1, #curvec do
    for k=1,3 do
      def_curr_wdg(k, curvec[i][k], 0)
    end

    calc_field_single({
        maxit=m.num_nonl_it, maxcop=m.error_perm, -- err_perm in %
        permode='restore'})
    dRR = 0
    RR = 0
%if model.get('sim_demagn', 0):
    repeat
       calc_field_single(1, 'actual', 0.005)
       --calc_field_single(300, 'reset', 0.005)
       stat, RR, dRR = calc_demag(5, m.hc_min)
       --printf("%g %g", RR, dRR)
    until math.abs(dRR) < 1e-5
%endif
    psi = {}
    for k=1,3 do
      psir, psii = flux_winding_wk(k)
      psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}
    end

    fr, ft, tq, fx, fy = force_torque()
    tqvec[i] = tq
    psivec[i] = psi
    rrvec[i] = RR
  end
  return psivec, tqvec, rrvec
end

%if type(model['curvec']) is list:
curamp = {${','.join([str(x) for x in model['curvec']])}} -- A
% else:
curamp = {${model['curvec']}} -- A
% endif
a=${model.get('num_par_wdgs', 1)}   -- parallel branches

curvec = {}
for i=1, #curamp do
    amp = math.sqrt(2)*curamp[i]/a
    curvec[i] = {amp, -amp, 0}
 end

ksym = m.num_poles/m.npols_gen

if num_agnodes ~= nil then
  dphi = 360/num_agnodes -- ndst[2] -- deg
else
  post_models("nodedistance", "ndst" )
  dphi = ndst[2] -- deg
end
nodes = math.floor(360/m.num_poles/dphi+0.5)
printf("Nodes in airgap total %g, Nodes per pole: %d", 360/dphi, nodes)
-- find a valid number of steps for a rotation:
nrot = nodes
while( nrot%2) == 0 do
  nrot = nrot//2
end
Q1 = get_dev_data("num_slots")
p = m.num_poles//2
dphi = 360//gcd(Q1, p)/nrot
print(string.format(" rotation steps: %d  current steps: %d\n", nrot, #curvec))

bags = {}

phi = 0
-- initialize rotate
rotate({
    airgap = m.fc_radius,    -- air gap radius
    region = "inside",       -- region to rotate
    mode   = "save"         -- save initial model state
})

file_psi = io.open("psi-torq-rem-rot.dat","w")
for n=1,nrot+1 do
  psi, tq = calc_flux_torq(phi, curvec)
  for i=1, #curvec do
    file_psi:write(string.format("%g ", phi))
    for k=1, 3 do
      file_psi:write(string.format("%g ", a*curvec[i][k]))
    end
    for k=1, 3 do
      file_psi:write(string.format("%g ", psi[i][k][1]))
    end
    file_psi:write(string.format("%g ", tq[i]))
    file_psi:write("\n")
  end

  phi = n*dphi
  rotate({angle=phi, mode="absolute"})
end
rotate({mode = "reset"})  -- restore the initial state (discard any changes)
file_psi:close()
