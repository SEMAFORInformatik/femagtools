-- calculate noload flux and airgap flux density
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
  f1=0
  a=${model.get('num_par_wdgs', 1)}   -- parallel branches

  for i=1, #curvec do
    bag[i] = {}
    amp = math.sqrt(2)*curvec[i]/a
    for k=1,3 do
    --Set currents
      alfa = (k-1)*2*math.pi/3
      def_curr_wdg(k, amp*math.cos(alfa), 0)
    end

  calc_field_single({
        maxit=300, maxcop=m.error_perm, -- err_perm in %
        permode='restore'})

  psi[i] = 0
  for k=1,3 do
    psi_re, psi_im = flux_winding_wk(k)
    alfa = (k-1)*2*math.pi/3
    psi[i] = psi[i] + ksym/a*m.arm_length*psi_re*math.cos(alfa)
  end

    post_models("induct(x)","b")    -- Calculate field distribution
      k = 1
      for j = 1, table.getn(b), 3 do
        pos[k] = b[j]
        bag[i][k] = b[j+1] -- radial component only
        k = k+1
      end
  end
  rec[1] = string.format("%g ", phi)
  file_psi:write(rec[1])
  for i=1, #curvec do
     rec[i+1] = string.format("%g %g ", curvec[i], math.sqrt(2)/3*psi[i])
     file_psi:write(rec[i+1])
  end
  print(table.concat(rec))
  file_psi:write("\n")
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
-- find a valid number of steps for a rotation of one rotor slot:
for n=5, nodes//2 do
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
file_psi = io.open("psi-rot-mag.dat","w")
for n=1,nrot+1 do

  pos, bag = calc_field(phi, curvec, file_psi)
  bags[n] = bag

  phi=dphi*n
  if (phi>=360.0/ksym) then
    rotate({airgap=m.fc_radius,
            angle=dphi-phi,
            region="inside",
            mode="increment"})  -- rotation pos outside model
  else
    rotate({airgap=m.fc_radius,
            angle=dphi,
            region="inside",
            mode="increment"})  -- rotation pos inside model
  end
end
rotate({
         mode = "reset"  -- restore the initial state (discard any changes)
       })
file_psi:close()

for i=1, #curvec do
  bagfile = string.format("noloadbag-%d.dat", i)
  file_bag = io.open(bagfile,"w")
  for k = 1, #pos do
     file_bag:write(string.format("%g ", pos[k]))
     for n=1, nrot+1 do
      file_bag:write(string.format("%g ",
         bags[n][i][k]))
     end
     file_bag:write("\n")
  end
  file_bag:close()
end
