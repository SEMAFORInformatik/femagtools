-- calculate noload flux and airgap flux density (magstatic mode)
--
-- model:
--    Q2 (number of rotor slots)
--    currvec (A rms)
--    num_par_wdgs (number of parallel winding groups, default 1)
--
function gcd(a,b)
     return b==0 and a or gcd(b,a%b)
end

function calc_field(phi, curvec, file_psi)
  psi={}
  pos={}
  bag={}
  rec={}
  cur={}
  f1=0
  alfa0 = math.pi/3
  a=${model.get('num_par_wdgs', 1)}   -- parallel branches

  for i=1, #curvec do
    bag[i] = {}
    amp = math.sqrt(2)*curvec[i]/a
    for k=1,3 do
    --Set currents
      alfa = (k-1)*2*math.pi/3 + alfa0
      def_curr_wdg(k, amp*math.cos(alfa), 0)
    end

    calc_field_single({
        maxit=300, maxcop=m.error_perm, -- err_perm in %
        permode='restore'})

    --psi[i] = 0
    for k=1,3 do
      psir, psii = flux_winding_wk(k)
      curr, curi = get_wdg_data("cur", k)
      cur[k] = {a*curr, a*curi}
      --alfa = (k-1)*2*math.pi/3
      --psi[i] = psi[i] + ksym/a*m.arm_length*psi_re*math.cos(alfa)
      psi[k] = {ksym*psir/a*m.arm_length, ksym*psii/a*m.arm_length}
    end

    post_models("induct(x)","b")    -- Calculate field distribution
    k = 1
    for j = 1, table.getn(b), 3 do
      pos[k] = b[j]
      bag[i][k] = b[j+1] -- radial component only
      k = k+1
    end

    file_psi:write(string.format("%g ", phi))
    for k=1, 3 do
      file_psi:write(string.format("%g %g ",
        cur[k][1], cur[k][2]))
    end
    for k=1, 3 do
      file_psi:write(string.format("%g %g ",
        psi[k][1], psi[k][2]))
    end
    file_psi:write("\n")
  end
  return pos, bag
end

ksym = m.num_poles/m.npols_gen
maxit = m.num_nonl_it
du_u0 = m.error_perm
rmod = m.perm_mode
cmod = 0

phi = 0.0
if num_agnodes ~= nil then
  dphi = 360/num_agnodes -- ndst[2] -- deg
else
  post_models("nodedistance", "ndst" )
  dphi = ndst[2] -- deg
end
nodes = math.floor(360/Q2/dphi+0.5)
printf("Nodes in airgap total %g, Nodes per rotor slot: %d", 360/dphi, nodes)
-- find a valid number of steps for a rotation of one rotor slot:
for n=nodes//4, nodes//2 do
  if nodes % n == 0 then
    dphi = nodes/n*dphi
    nrot = n
    break
  end
end
if nrot == nil then
  nrot = nodes
end
curvec = {${','.join([str(x) for x in model['curvec']])}} -- A rms

print("\nNo load flux simulation (DC) with rotation\n")
print(string.format(" rotation steps: %d  current steps: %d\n", nrot, #curvec))

bags = {}

phi = 0
-- initialize rotate
rotate({
    airgap = m.fc_radius,    -- air gap radius
    region = "inside",       -- region to rotate
    mode   = "save"         -- save initial model state
})

file_psi = io.open("psi-rot-mag.dat","w")
for n=1,nrot+1 do

  pos, bag = calc_field(phi, curvec, file_psi)
  bags[n] = bag
  phi = n*dphi
  rotate({angle=phi, mode="absolute"})
end
rotate({mode = "reset"})  -- restore the initial state (discard any changes)
file_psi:close()

for i=1, #curvec do
  bagfile = string.format("noloadbag-%d.dat", i)
  file_bag = io.open(bagfile,"w")
  for k = 1, #pos do
     file_bag:write(string.format("%g ", pos[k]))
     for n=1, nrot+1 do
      file_bag:write(string.format("%g ",
         bags[n][i][k]))  -- Br, rotpos, cur, pos
     end
     file_bag:write("\n")
  end
  file_bag:close()
end
