
-- Switched reluctance rotor

r = ${model['r']}*1e3    -- rounding radius
alfa = ${model['alfa']}  -- pole width factor 0 < alfa < 1

P = m.num_poles
Pm = m.npols_gen
m.zeroangl = 0  -- zero angle in degree
taup = m.zeroangl + 360/P -- pole pitch in degree


fml = require ("fml")

h = (da2 - dy2)/2
M = fml.Point:Create(0, 0)
C1 = fml.Circle:Create(M, dy2/2)
C1x = fml.Circle:Create(M, da2/2-h/4)
L1 = fml.Line:Create(M, 0)
P1 = fml.Point:Intersection(L1, C1, 2)
C2 = fml.Circle:Create(M, da2/2)
P2 = fml.Point:Intersection(L1, C2, 2)
L2 = fml.Line:Create(M, alfa*taup/2)
L3 = fml.Line:Create(M, taup/2)
P3 = fml.Point:Intersection(L2, C2, 2)
L4 = fml.Line:Parallel(L1, P3)
P4 = fml.Point:Intersection(L3, C2, 2)
P5 = fml.Point:Intersection(L4, L3)
P6 = fml.Point:Intersection(L3, C1, 2)
L5 = fml.Line:Create(P5, taup)

M1, P7, P8 = fml.Point:Rounding(L4, L5, r, 1)
C3 = fml.Circle:Create(M1, r)
P9 = fml.Point:Intersection(L3, C3, 1)

P3x = fml.Point:Intersection(L4, C1x, 2)

ndt(agndst)
nc_circle(P2.x, P2.y, P3.x, P3.y, 0)
nc_circle(P3.x, P3.y, P4.x, P4.y, 0)
nc_line(P3.x, P3.y, P3x.x, P3x.y, 0)

ndt(4*agndst)

nc_line(P1.x, P1.y, P2.x, P2.y, 0)
nc_line(P4.x, P4.y, P9.x, P9.y, 0)
nc_line(P3x.x, P3x.y, P8.x, P8.y, 0)
nc_circle_m(P9.x, P9.y, P8.x, P8.y, M1.x, M1.y, 0)
nc_line(P9.x, P9.y, P6.x, P6.y, 0)
nc_circle(P1.x, P1.y, P6.x, P6.y, 0)
create_mesh_se((P3.x + P9.x)/2, (P2.y + P3.y)/2)
def_new_subreg((P3.x + P9.x)/2, (P2.y + P3.y)/2, "Rotor", "skyblue")
create_mesh_se(P4.x, (P3.y + P4.y)/2)
 
mirror_nodechains(P4.x, P4.y, P6.x, P6.y)
x1, y1 = pd2c(da2/2, taup)
x2, y2 = pd2c(dy2/2, taup)

rotate_copy_nodechains(P1.x, P1.y, P2.x, P2.y, x1, y1, x2, y2, m.npols_gen-1)

if mcvkey_yoke ~= nil and mcvkey_yoke ~= 'dummy' then
  m.rlength         =     ${model.get('rlength', 1)*100}  
  def_mat_fm_nlin((P3.x + P5.x)/2, (P2.y + P3.y)/2, "blue", mcvkey_yoke, m.rlength)
else
  def_mat_fm((P3.x + P5.x)/2, (P2.y + P3.y)/2, 1000, m.rlength)
end

if(m.npols_gen == m.num_poles) then
  def_bcond_vpo(P1.x, P1.y, -P1.x, P1.y)
  def_bcond_vpo(-P1.x, P1.y, P1.x, P1.y)
end
