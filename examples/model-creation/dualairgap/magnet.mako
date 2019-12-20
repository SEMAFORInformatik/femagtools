
hm = ${model['hm']}*1e3
alpha_m = ${model['alpha_m']}

x = {}
y = {} 
for i=1, 12 do
  x[i]=0
  y[i]=0
end

hrs = (dy2/2-da1[2]/2)-hm

x[1], y[1] = dy2/2, 0                                  -- Yoke outside
x[2], y[2] = pd2c(dy2/2, 180/m.num_poles)
x[3], y[3] = pd2c(dy2/2-hrs, 180/m.num_poles)
x[4], y[4] = dy2/2-hrs, 0

x[5], y[5] = pd2c(dy2/2-hrs, 180/m.num_poles*alpha_m)   -- Permanent magnet
x[6], y[6] = pd2c(da1[2]/2, 180/m.num_poles*alpha_m)
x[7], y[7] = da1[2]/2, 0 
x[8], y[8] = pd2c(da1[2]/2, 360/(m.num_poles*2))

x[9], y[9] = da1[1]/2, 0                             -- Yoke inside
x[10], y[10] = pd2c(da1[1]/2, 180/m.num_poles)
x[11], y[11] = pd2c(dy1/2, 180/m.num_poles)
x[12], y[12] = dy1/2, 0

ndt(1)

nc_circle(x[1],y[1],x[2],y[2],0)             -- Yoke outside
nc_line(x[2],y[2],x[3],y[3],0)
nc_circle(x[4],y[4],x[5],y[5],0)
if alpha_m ~= 1.0 then
  nc_circle(x[5],y[5],x[3],y[3],0)
end
nc_line(x[4],y[4],x[1],y[1],0)
x0, y0 = pd2c(dy2/2-hrs/2, 90/m.num_poles)
create_mesh_se(x0, y0)
def_new_subreg(x0, y0, "Iron", blue)

nc_line(x[5],y[5],x[6],y[6],0)               -- Permanent magnet
nc_circle(x[7],y[7],x[6],y[6],n2*alpha_m)
nc_line(x[7],y[7],x[4],y[4],0)
if alpha_m ~= 1.0 then
  nc_circle(x[6],y[6],x[8],y[8],n2*(1-alpha_m)+1)
  nc_line(x[8],y[8],x[3],y[3],0)
end
x0, y0 = pd2c(da1[2]/2+hm/2, (180/m.num_poles*alpha_m)/2)
create_mesh_se(x0, y0)

x0, y0 = pd2c(da1[2]/2+hm/2, (180/m.num_poles*(1+alpha_m)/2))
create_mesh_se(x0, y0)


nc_circle(x[9],y[9],x[10],y[10],n1)         -- Yoke inside
nc_line(x[10],y[10],x[11],y[11],0)
nc_circle(x[12],y[12],x[11],y[11],0)
nc_line(x[12],y[12],x[9],y[9],0)
x0, y0 = pd2c((da1[1]+dy1)/4, 90/m.num_poles)
create_mesh_se(x0, y0)
add_to_subreg(x0, y0, "Iron")

x[13], y[13] = da1[1]/2+1/3*ag[1], 0   -- airgap 1
x[14], y[14] = pd2c(x[13], 180/m.num_poles)
nc_circle(x[13],y[13],x[14],y[14],n1)
nc_line(x[9], y[9], x[13], y[13], 0)
nc_line(x[10], y[10], x[14], y[14], 0)

x0, y0 = pd2c( da1[1]/2 + 1/6*ag[1], 180/m.num_poles)
create_mesh_se(x0, y0)
mirror_nodechains(x[14],y[14],x[11],y[11]) -- Kopieren vom 1/2 Teilsegment

x[15], y[15] = da1[2]/2-1/3*ag[2], 0   -- airgap 2
x[16], y[16] = pd2c(x[15], 180/m.num_poles)
nc_circle(x[15],y[15],x[16],y[16],n2)
nc_line(x[8], y[8], x[16], y[16], 0)
nc_line(x[7], y[7], x[15], y[15], 0)

x0, y0 = pd2c( da1[2]/2 - 1/6*ag[1], 90/m.num_poles)
create_mesh_se(x0, y0)

mirror_nodechains(x[2],y[2],x[16],y[16]) -- Kopieren vom 1/2 Teilsegment

-- Material properties
relperm = 1000
x0, y0 = pd2c(dy2/2-hrs/2, 90/m.num_poles)
def_mat_fm(x0, y0, yellow, relperm, 100)      -- Yoke outside
xm, ym = pd2c(da1[2]/2 + hm/2, 90/m.num_poles)
def_mat_pm(xm, ym, magenta, m.remanenc, m.relperm, 0, m.radial, 100)     -- Permanentmagnet 1
xm, ym = pd2c(da1[2]/2 + hm/2, 270/m.num_poles)
def_mat_pm(xm, ym, cyan, m.remanenc, m.relperm, 180, m.radial, 100)      -- Permanentmagnet 2
x0, y0 = pd2c((da1[1]+dy1)/4, 90/m.num_poles)
def_mat_fm(x0, y0, yellow, relperm, 100)      -- Rueckschluss innen
