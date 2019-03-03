
s = ${model['slot_width']*1e3}
bz = ${model['tooth_width']*1e3}
h1 = ${model['slot_h1']*1e3}
h2 = ${model['slot_h2']*1e3}
hj = ${model['yoke_height']*1e3}

x = {}
y = {} 
for i=1, 15 do
  x[i]=0
  y[i]=0
end

Da = dy1
Di = da1
Q = m.tot_num_slot

x[1],y[1] = pd2c(Da/2,0)
x[2],y[2] = pd2c(Da/2,180/Q)
x[3] = Di/2*math.cos(math.asin(s/Di))
y[3] = s/2
x[4] = x[3]+h1
y[4] = s/2
x[5],y[5] = pd2c(Di/2,180/Q)
x[6] = Di/2+h1+h2;
y[6] = y[5]/x[5]*x[6]-bz/2/math.cos(pi/Di)
x[7] = Da/2-hj;
y[7] = y[5]/x[5]*x[7]-bz/2/math.cos(pi/Di)
x[8] = x[7]
x[9],y[9] = pd2c(vlen(x[4],y[4]),180/Q)
x[10] = (y[6]+x[5]/y[5]*x[6])/(y[5]/x[5]+x[5]/y[5])
y[10] = y[5]/x[5]*x[10]
x[11] = (y[7]+x[5]/y[5]*x[7])/(y[5]/x[5]+x[5]/y[5])
y[11] = y[5]/x[5]*x[11]
x[12] = Di/2
x[13] = x[4]
x[14] = Di/2-ag/3
x[15],y[15] = pd2c(Di/2-ag/3,180/Q)

-- create node chains

agnp = m.nodedist
agndst = ag*2/3
ndt(agndst)

nc_circle(x[14],y[14],x[15],y[15],360/Q/2/agnp+1) 
nc_circle(x[1],y[1],x[2],y[2],0)
nc_circle(x[13],y[13],x[4],y[4],0)
nc_circle(x[3],y[3],x[5],y[5],0)
nc_line(x[3],y[3],x[4],y[4],0)
nc_line_cont(x[6],y[6],0)
nc_line_cont(x[7],y[7],0)
nc_line_cont(x[8],y[8],0)
nc_line(x[12],y[12],x[13],y[13],0)
nc_line_cont(x[8],y[8],0)
nc_line_cont(x[1],y[1],0)
nc_line(x[14],y[14],x[12],y[12],0)
nc_line(x[15],y[15],x[5],y[5],0)
nc_line_cont(x[9],y[9],0)
nc_line_cont(x[10],y[10],0)
nc_line_cont(x[11],y[11],0)
nc_line_cont(x[2],y[2],0)

-- Meshing

mesh.con1 = 0.1                   -- Vernetzungssteuerung

create_mesh_se(Da/2-hj/2,0+hj/2)
create_mesh_se((Da+Di)/4,s/4)
create_mesh_se(Di/2+h1/2,s/4)

-- Definition of Subregions

def_new_subreg(Da/2-hj/2,0+hj/2,"Stat",blue)

-- Mirror and Rotate

mirror_nodechains(x[2],y[2],x[15],y[15])

x0,y0 = pd2c(Di/2-ag/3,0)
x1,y1 = pd2c(Da/2,0)
x2,y2 = pd2c(Da/2,360/Q)
x3,y3 = pd2c(Di/2-ag/3,360/Q)
rotate_copy_nodechains(x0,y0,x1,y1,x2,y2,x3,y3,m.num_slots-1)

-- winding location

tauq = 360/Q               -- Nutteilungswinkel
Rq = (Di/2+Da/2-hj)/2      -- mittlerer Nutradius
m.wdg_location = 1
m.xcoil_1, m.ycoil_1 = pd2c(Rq,tauq/4)
m.xcoil_2, m.ycoil_2 = pd2c(Rq,tauq-tauq/4)

-- Material properties
rlen=100
if mcvkey_yoke ~= 'dummy' then
  def_mat_fm_nlin(Da/2-hj/2, 1.0, blue, mcvkey_yoke, rlen)
else
  urr = 1000
  def_mat_fm(Da/2-hj/2,ag,urr,rlen)
end
