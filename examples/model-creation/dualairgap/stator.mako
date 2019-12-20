beta_s = ${model['beta_s']}
nodepitch1 = ${model['nodepitch1']}
nodepitch2 = ${model['nodepitch2']} 

x = {}
y = {}
for i=1, 6 do
  x[i]=0
  y[i]=0
end

x[1], y[1] = da1[2]/2-ag[2], 0    
x[2], y[2] = pd2c(da1[2]/2-ag[2], (180/m.num_poles)*(1-beta_s))
x[3], y[3] = pd2c(da1[1]/2+ag[1], (180/m.num_poles)*(1-beta_s))
x[4], y[4] = da1[1]/2+ag[1], 0
x[5], y[5] = pd2c(da1[2]/2-ag[2], 180/m.num_poles)
x[6], y[6] = pd2c(da1[1]/2+ag[1], 180/m.num_poles)

n1 = 180/m.num_poles/nodepitch1+1  -- Anzahl Knoten von Luftspalt 1
n1_S1 = n1*(1-beta_s)              -- Anzahl Knoten von Luftspalt 1 an linker Seite von Spule
n1_S2 = n1-n1_S1+1                           -- Anzahl Knoten von Luftspalt 1 oberhalb der lin

n2 = 180/m.num_poles/nodepitch2+1  -- Anzahl Knoten von Luftspalt 2
n2_S1 = n2*(1-beta_s)              -- Anzahl Knoten von Luftspalt 2 an rechter Seite von Spule
n2_S2 = n2-n2_S1+1                 -- Anzahl Knoten von Luftspalt 2 oberhalb der rechten Seite der Spule


nc_circle(x[1],y[1],x[2],y[2],n2_S1)           -- Winding area
nc_line(x[2],y[2],x[3],y[3],0)
if beta_s ~= 0.0 then
  nc_circle(x[2],y[2],x[5],y[5],n2_S2)
  nc_circle(x[3],y[3],x[6],y[6],n1_S2)
  nc_line(x[6],y[6],x[5],y[5],0)
end
nc_circle(x[4],y[4],x[3],y[3],n1_S1)
nc_line(x[4],y[4],x[1],y[1],0)

x0, y0 = pd2c((da1[1] + da1[2])/4, (180/m.num_poles)*(1-beta_s)/2)
create_mesh_se(x0, y0)
def_new_subreg(x0,y0,"YK",blue)
if beta_s ~= 0.0 then
  x0, y0 = pd2c((da1[1] + da1[2])/4, 180/m.num_poles*(1-beta_s/2))
  m.xcoil_1, m.ycoil_1 = x0, y0
  create_mesh_se(x0, y0)
  x0, y0 = pd2c((da1[1] + da1[2])/4, 180/m.num_poles*(1+beta_s/2))
  m.xcoil_2, m.ycoil_2 = x0, y0
end

x[7], y[7] = da1[1]/2+2/3*ag[1], 0          -- inner airgap 
x[8], y[8] = pd2c(x[7], 180/m.num_poles)
nc_circle(x[7],y[7],x[8],y[8],n1)
nc_line(x[4], y[4], x[7], y[7], 0)
nc_line(x[6], y[6], x[8], y[8], 0)

x0, y0 = pd2c( da1[1]/2 + 5/6*ag[1], 180/m.num_poles)
create_mesh_se(x0, y0)

x[9], y[9] = da1[2]/2-2/3*ag[2], 0           -- outer airgap
x[10], y[10] = pd2c(x[9], 180/m.num_poles)
nc_circle(x[9],y[9],x[10],y[10],n2)
nc_line(x[1], y[1], x[9], y[9], 0)
nc_line(x[5], y[5], x[10], y[10], 0)

x0, y0 = pd2c( da1[2]/2 - 5/6*ag[2], 90/m.num_poles)
create_mesh_se(x0, y0)

mirror_nodechains(x[10],y[10],x[8],y[8]) -- Copy segment

-- Material properties
x0, y0 = pd2c((da1[1] + da1[2])/4, (180/m.num_poles)*(1-beta_s)/2)
def_mat_fm(x0, y0, 1.0001, 100)
