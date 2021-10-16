
magn_height = ${model['magn_height']*1e3}    -- magnet height [m]
rel_magn_width = ${model['rel_magn_width']}    -- rel magnet width
magn_length = ${model['magn_length']*1e3}    -- magnet length [m]
mag_ori = '${model['magn_ori']}'
bridge_height = ${model['bridge_height']*1e3}


dair_hull = 1.2*dy1

yoke_height = dy1/2-da1/2-magn_height

rlenPM = 100*(magn_length/m.arm_length) -- length of PM in %

m.zeroangl = 0

fml = require("fml")

-- Berechnung der Koordinaten
M = fml.Point:Create(0,0)
C1 = fml.Circle:Create(M, da1/2)
C2 = fml.Circle:Create(M, da1/2+bridge_height)
C3 = fml.Circle:Create(M, da1/2+magn_height)
C4 = fml.Circle:Create(M, dy1/2)
C5 = fml.Circle:Create(M, dair_hull/2)
taup = 360/m.num_poles
L1 = fml.Line:Create(M, m.zeroangl)
L2 = fml.Line:Create(M, m.zeroangl+taup/2)
L3 = fml.Line:Create(M, m.zeroangl+(1-rel_magn_width)*taup/2)

P1 = fml.Point:Intersection(L1, C1, 2)
P2 = fml.Point:Intersection(L1, C2, 2)
P3 = fml.Point:Intersection(L1, C3, 2)
P4 = fml.Point:Intersection(L1, C4, 2)
P5 = fml.Point:Intersection(L1, C5, 2)
P6 = fml.Point:Intersection(L2, C5, 2)
P7 = fml.Point:Intersection(L2, C4, 2)
P8 = fml.Point:Intersection(L2, C3, 2)
P9 = fml.Point:Intersection(L2, C2, 2)
P10 = fml.Point:Intersection(L2, C1, 2)
P11 = fml.Point:Intersection(L3, C1, 2)
P12 = fml.Point:Intersection(L3, C2, 2)
P13 = fml.Point:Intersection(L3, C3, 2)

-- Knotenketten Staender
ndt(agndst)
nc_circle(P1.x, P1.y, P11.x, P11.y, 0)
nc_circle(P11.x, P11.y, P10.x, P10.y, 0)
nc_circle(P2.x, P2.y, P12.x, P12.y, 0)
nc_circle(P3.x, P3.y, P13.x, P13.y, 0)
nc_circle(P13.x, P13.y, P8.x, P8.y, 0)

nc_line(P1.x, P1.y, P2.x, P2.y,0)
nc_line_cont(P3.x, P3.y, 0)
nc_line(P11.x, P11.y, P12.x, P12.y, 0)
nc_line_cont(P13.x, P13.y, 0)
nc_line(P10.x, P10.y, P9.x, P9.y, 0)
nc_line_cont(P8.x, P8.y, 0)

ndt(2*agndst)
nc_circle(P4.x, P4.y, P7.x, P7.y, 0)
nc_line(P3.x, P3.y, P4.x, P4.y, 0)
nc_line(P8.x, P8.y, P7.x, P7.y, 0)
ndt(3*agndst)
nc_line(P4.x, P4.y, P5.x, P5.y, 0)
nc_line(P7.x, P7.y, P6.x, P6.y, 0)
ndt(4*agndst)
nc_circle(P5.x, P5.y, P6.x, P6.y, 0)

-- Meshing 

-- Air outside
x, y = pd2c(dy1/2 + 0.1, m.zeroangl + taup/4)
create_mesh_se(x, y)
-- Air inside
x, y = pd2c(da1/2 + 0.1, m.zeroangl + taup/4*(1-rel_magn_width))
create_mesh_se(x, y)

-- Yoke
x, y = pd2c(dy1/2 - yoke_height/2, m.zeroangl + taup/4)
create_mesh_se(x, y)
def_new_subreg(x, y, "styk",blue)
x, y = pd2c(dy1/2 - yoke_height - bridge_height/2, m.zeroangl + taup/4*(1-rel_magn_width))
create_mesh_se(x, y)
add_to_subreg(x, y, "styk",blue)

-- PM
x, y = pd2c(da1/2+magn_height/2,taup/4)
create_mesh_se(x, y)

mirror_nodechains(P6.x, P6.y, P10.x, P10.y)
x1, y1 = pd2c(dair_hull/2, m.zeroangl+taup)
x2, y2 = pd2c(da1/2, m.zeroangl+taup)
rotate_copy_nodechains(P1.x,P1.y,P5.x,P5.y,
                       x1,y1,x2,y2,m.npols_gen-1)

------------------------------
----- materials definition ---------------
------------------------------
if mcvkey_yoke ~= 'dummy' then
  def_mat_fm_nlin(dy1/2-yoke_height/2, ag/3, "blue", mcvkey_yoke, rlenPM)
else
  def_mat_fm(dy1/2-yoke_height/2, ag/3, "blue", 1000, rlenPM)
end
-- Permanentmagnete
for i=0, m.npols_gen-1, 2 do
	alfa = i*taup + (1-rel_magn_width/2)*taup/2
	x,y = pd2c(da1/2+magn_height/2,alfa)
	if mag_ori=="radial" then
        gam = 0
	else
        gam = alfa
	end
	def_mat_pm(x,y,"red",m.remanenc,m.relperm,gam,mag_ori,rlenPM)
	alfa = alfa + rel_magn_width*taup/2
	x,y = pd2c(da1/2+magn_height/2,alfa)
	if mag_ori=="radial" then
        gam = 0
	else
        gam = alfa
	end
	def_mat_pm(x,y,"red",m.remanenc,m.relperm,gam,mag_ori,rlenPM)
end

for i=1, m.npols_gen-1, 2 do
	alfa = i*taup + (1-rel_magn_width/2)*taup/2
	x,y = pd2c(da1/2+magn_height/2,alfa)
	if mag_ori=="radial" then
        gam = 180
	else
        gam = alfa-180
	end
	def_mat_pm(x,y,"green",m.remanenc,m.relperm,gam,mag_ori,rlenPM)
	alfa = alfa + rel_magn_width*taup/2
	x,y = pd2c(da1/2+magn_height/2,alfa)
	if mag_ori=="radial" then
        gam = 180
	else
        gam = alfa-180
	end
	def_mat_pm(x,y,"green",m.remanenc,m.relperm,gam,mag_ori,rlenPM)
end
