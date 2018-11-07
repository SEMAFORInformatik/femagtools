
-- Model parameters

ds = ${model['shaft_diam']}*1e3
hm = ${model['magn_height']}*1e3
bm = ${model['magn_width']}*1e3
ws = ${model['slot_width']}*1e3

P = m.num_poles
Pm = m.npols_gen
m.zeroangl = 0  -- zero angle in degree
taup = m.zeroangl + 2*360/P -- pole pitch in degree

---------
-- FML --
---------
fml = require("fml")

fml.SetLocalCS(0, 0, 0)

P = {} -- Points
L = {} -- Lines
C = {} -- Circles

M = fml.Point:Create(0, 0)

C[1] = fml.Circle:Create(M, ds/2)
C[2] = fml.Circle:Create(M, dy2/2)
C[3] = fml.Circle:Create(M, dy2/2+bm)
C[4] = fml.Circle:Create(M, da2/2)

L[1] = fml.Line:Create(M, m.zeroangl)
L[2] = fml.Line:Create(M, m.zeroangl + taup/4)
L[3] = fml.Line:Parallel(L[2], hm/2)
L[4] = fml.Line:Parallel(L[2], ws/2)

P[8] = fml.Point:Intersection(L[2], C[2], 2)
P[9] = fml.Point:Intersection(L[2], C[3], 2)
L[5] = fml.Line:Create(P[9], m.zeroangl + taup/4 + 90)
L[6] = fml.Line:Create(P[8], m.zeroangl + taup/4 + 90)

P[3] = fml.Point:Intersection(L[1], C[4], 2)
P[4] = fml.Point:Intersection(L[4], C[4], 2)
P[5] = fml.Point:Intersection(L[4], L[5])
P[6] = fml.Point:Intersection(L[5], L[3])
P[7] = fml.Point:Intersection(L[6], L[3])
P[10] = fml.Point:Intersection(L[2], C[4], 2)

L[7] = fml.Line:Create(P[7], m.zeroangl + 90)
P[1] = fml.Point:Intersection(L[1], L[7])

r, phi = c2pd(P[6].x, P[6].y)
C[5] = fml.Circle:Create(M, r)
P[2] = fml.Point:Intersection(L[1], C[5], 2)

-----------------
-- Node chains --
-----------------
ndt(agndst)
nc_circle(P[3].x, P[3].y, P[4].x, P[4].y, 0)
nc_circle(P[4].x, P[4].y, P[10].x, P[10].y, 0)

ndt(1.4*agndst)
nc_line(P[3].x, P[3].y, P[2].x, P[2].y, 0)
nc_line(P[10].x, P[10].y, P[9].x, P[9].y, 0)
nc_line_cont(P[5].x, P[5].y, 0)
nc_line(P[4].x, P[4].y, P[5].x, P[5].y, 0)
nc_line_cont(P[6].x, P[6].y, 0)

ndt(2*agndst)
nc_circle(P[2].x, P[2].y, P[6].x, P[6].y, 0)
nc_line(P[2].x, P[2].y, P[1].x, P[1].y, 0)
nc_line_cont(P[7].x, P[7].y, 0)
nc_line(P[6].x, P[6].y, P[7].x, P[7].y, 0)
nc_line_cont(P[8].x, P[8].y, 0)
nc_line_cont(P[9].x, P[9].y, 0)

x0, y0 = pd2c((da2/2 + dy2/2+bm)/2, taup/4 + m.zeroangl - 0.1)
create_mesh_se(x0, y0)

r, phi = c2pd((P[2].x + P[3].x)/2, (P[2].y + P[3].y)/2)
x0, y0 = pd2c(r, m.zeroangl + 0.1)
create_mesh_se(x0, y0)
def_new_sreg(x0, y0, "RotorYoke", yellow)
x0, y0 = (P[1].x+P[2].x)/2, (P[1].y+P[7].y)/2
create_mesh_se(x0, y0)
add_to_sreg(x0, y0,"RotorYoke")

x0, y0 = pd2c(dy2/2+bm/2, taup/4 + m.zeroangl -0.1)
create_mesh_se(x0, y0)
def_new_sreg(x0, y0, "Magnet", red)

mirror_nodechains(P[10].x, P[10].y, P[8].x, P[8].y)

x2, y2 = pd2c(da2/2, taup/2 + m.zeroangl)
r, phi = c2pd(P[1].x, P[1].y)
x3, y3 = pd2c(r, taup/2 + m.zeroangl)
rotate_copy_nodechains(P[1].x, P[1].y, P[3].x, P[3].y, x2, y2, x3, y3, Pm-1)

-- -----------------------------
-- ---- Material Properties ----
-- -----------------------------

-- --  PM yoke
urr=1000
x0, y0 = (da2+dy2)/4, 0.1
if mcvkey_yoke ~= 'dummy' then
   def_mat_fm_nlin(x0, y0, blue, mcvkey_yoke, 100)
else
   def_mat_fm(x0, y0, 1000.0, 100)
end

-- -- Permanent magnets
for i=0, m.npols_gen-1 do
   phi = (2*i+1)*(taup/4 + m.zeroangl)
   x0,y0 = pr2c((da2+dy2)/4, phi*math.pi/180 - math.atan2(hm/4,(da2+dy2)/4))
   x1,y1 = pr2c((da2+dy2)/4, phi*math.pi/180 + math.atan2(hm/4,(da2+dy2)/4))
   if ( i % 2 == 0 ) then
      def_mat_pm(x0, y0,red,m.remanenc, m.relperm, phi+90, m.parallel,100)
      def_mat_pm(x1, y1,red,m.remanenc, m.relperm, phi+90, m.parallel,100)
   else
      def_mat_pm(x0, y0, green, m.remanenc, m.relperm, phi-90, m.parallel,100)
      def_mat_pm(x1, y1, green, m.remanenc, m.relperm, phi-90, m.parallel,100)
   end
end
