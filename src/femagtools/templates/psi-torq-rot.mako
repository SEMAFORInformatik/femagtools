-- calculate flux linkages and torque (magstatic mode)
--
-- model:
--    curvec (A) current samples (amplitudes)
--    num_par_wdgs (number of parallel winding groups, default 1)
--
-- creates file psi-torq-rot.dat in current directory with columns:
--  displ curr1 curr2 curr3 psi1 psi2 psi3 torq
--
function gcd(a, b)
	return b==0 and a or gcd(b,a%b)
end

function calc_flux_torq(phi, curvec)
  psivec={}
  tqvec={}
  for i=1, #curvec do
    for k=1,3 do
      def_curr_wdg(k, curvec[i][k], 0)
    end

    calc_field_single({
        maxit=m.num_nonl_it, maxcop=m.error_perm, -- err_perm in %
        permode='restore'})
    psi = {}
    for k=1,3 do
      psir, psii = flux_winding_wk(k)
      psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}
    end

    fr, ft, tq, fx, fy = force_torque()
    tqvec[i] = tq
    psivec[i] = psi
  end
  return psivec, tqvec
end

%if type(model['curvec']) is list:
curamp = {${','.join([str(x) for x in model['curvec']])}} -- A
% else:
curamp = {${model['curvec']}} -- A
% endif
a=${model.get('num_par_wdgs', 1)}   -- parallel branches

curvec = {}
for i=1, #curamp do
    amp = curamp[i]/a
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

phi = 0
-- initialize rotate
rotate({
    airgap = m.fc_radius,    -- air gap radius
    region = "inside",       -- region to rotate
    mode   = "save"         -- save initial model state
})

file_psi = io.open("psi-torq-rot.dat","w")
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
