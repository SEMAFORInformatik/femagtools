-- Model parameters

giron = ${model['gap_ma_iron']}*1e3 
hm = ${model['magn_height']}*1e3
bm = ${model['magn_width']}*1e3
gair = ${model['gap_air_iron']}*1e3
gam = ${model['gamma']}

r = da2/2 - giron - hm

P = m.num_poles
Pm = m.npols_gen

taup = 2*360/P -- pole pitch in degree

---------
-- FML --
---------
fml = require("fml")

fml.SetLocalCS(0, 0, 0)

P = {} -- Points
L = {} -- Lines
C = {} -- Circles

M = fml.Point:Create(0, 0)

L[1] = fml.Line:Create(M, m.zeroangl)
L[2] = fml.Line:Create(M, m.zeroangl + taup/4)

C[1] = fml.Circle:Create(M, dy2/2)
P[1] = fml.Point:Intersection(L[1], C[1], 2)
P[8] = fml.Point:Intersection(L[2], C[1], 2)

C[2] = fml.Circle:Create(M, r)
C[3] = fml.Circle:Create(M, da2/2 - gair)
P[5] = fml.Point:Intersection(L[2], C[3], 2)

C[4] = fml.Circle:Create(M, da2/2)
P[3] = fml.Point:Intersection(L[1], C[4], 2)
P[4] = fml.Point:Intersection(L[2], C[4], 2)

P[7] = fml.Point:Intersection(L[2], C[2], 2)
L[3] = fml.Line:Create(P[7], m.zeroangl + taup/4 + 90)
L[4] = fml.Line:Parallel(L[3], hm)
P[6] = fml.Point:Intersection(L[2], L[4])

L[5] = fml.Line:Parallel(L[2], bm/2)
P[10] = fml.Point:Intersection(L[4], L[5])
L[7] = fml.Line:Create(P[10], m.zeroangl + taup/4 - gam)
L[6] = fml.Line:Parallel(L[7], hm)
P[9] = fml.Point:Intersection(L[3], L[6])

m1, P[11], P[12] = fml.Point:Rounding(L[6], C[3], hm/2, 6)

m2, p1, P[13] = fml.Point:Rounding(L[7], C[3], hm/2, 7)

L[8] = fml.Line:Perpendicular(P[11], L[1])
P[2] = fml.Point:Intersection(L[1], L[8])
P[14] = fml.Point:Intersection(L[5], L[3])

-----------------
-- Node chains --
-----------------
ndt(agndst)
nc_circle(P[3].x, P[3].y, P[4].x, P[4].y, 0)

ndt(2*agndst)
nc_line(P[11].x, P[11].y, P[2].x, P[2].y, 0)
nc_line_cont(P[3].x, P[3].y, 0)
nc_line(P[7].x, P[7].y, P[6].x, P[6].y, 0)
nc_line_cont(P[4].x, P[4].y, 0)

ndt(3*agndst)
nc_line(P[7].x, P[7].y, P[9].x, P[9].y, 0)
nc_line_cont(P[11].x, P[11].y, 0)

nc_circle_m(P[11].x, P[11].y, P[12].x, P[12].y, m1.x, m1.y, 5)
nc_circle_m(P[12].x, P[12].y, P[13].x, P[13].y, m2.x, m2.y, 3)

nc_line(P[13].x, P[13].y, P[10].x, P[10].y, 0)
nc_line_cont(P[6].x, P[6].y, 0)
nc_line(P[14].x, P[14].y, P[10].x, P[10].y, 0)

ndt(4*agndst)
nc_line(P[1].x, P[1].y, P[2].x, P[2].y, 0)
nc_line(P[8].x, P[8].y, P[7].x, P[7].y, 0)

ndt(5*agndst)
nc_circle(P[1].x, P[1].y, P[8].x, P[8].y, 0)

-----
-- Meshing

x0, y0 = pd2c(da2/2 - gair/2, m.zeroangl + taup/8)
create_mesh_se(x0, y0)
def_new_sreg(x0, y0, "RotorYoke", mediumblue)

x0, y0 = pd2c((dy2/2 + r)/2, m.zeroangl + taup/8)
create_mesh_se(x0, y0)
add_to_sreg(x0, y0,"RotorYoke")

x0, y0 = pd2c(r + hm/2, m.zeroangl + taup/4 -0.1)
create_mesh_se(x0, y0)

create_mesh_se(m1.x, m1.y)

mirror_nodechains(P[4].x, P[4].y, P[8].x, P[8].y)

x2, y2 = pd2c(da2/2, taup/2 + m.zeroangl)
r1, phi = c2pd(P[1].x, P[1].y)
x3, y3 = pd2c(r1, taup/2 + m.zeroangl)
rotate_copy_nodechains(P[1].x, P[1].y, P[3].x, P[3].y, x2, y2, x3, y3, Pm-1)

-- -----------------------------
-- ---- Material Properties ----
-- -----------------------------

-- --  PM yoke
urr=1000
x0, y0 = pd2c((da2+dy2)/4, 0.1+m.zeroangl)
if mcvkey_yoke ~= 'dummy' then
   def_mat_fm_nlin(x0, y0, blue, mcvkey_yoke, 100)
else
   def_mat_fm(x0, y0, mediumblue, 1000.0, 100)
end

-- -- Permanent magnets
for i=0, m.npols_gen-1 do
   phi = (2*i+1)*taup/4 + m.zeroangl
   x0,y0 = pd2c(r + hm/2, phi+0.1)
   x1,y1 = pd2c(r + hm/2, phi-01)
   if ( i % 2 == 0 ) then
      def_mat_pm(x0, y0,red,m.remanenc, m.relperm, phi, m.parallel,100)
      def_mat_pm(x1, y1,red,m.remanenc, m.relperm, phi, m.parallel,100)
   else
      def_mat_pm(x0, y0, green, m.remanenc, m.relperm, phi+180, m.parallel,100)
      def_mat_pm(x1, y1, green, m.remanenc, m.relperm, phi+180, m.parallel,100)
   end
end
