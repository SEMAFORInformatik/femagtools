
------------------------------------------------------------------
-- Stator --------------------------------------------------------
------------------------------------------------------------------

Q = 12          
P = 10          

Da = dy1        
Di = da1         

s = ${model['sw']*1e3}           -- slot opening
bz = ${model['tw']*1e3}          -- tooth width
h1 = ${model['slot_h1']*1e3}        -- tooth height 1
h2 = ${model['slot_h2']*1e3}          -- tooth height 2
hj = 8          -- yoke height
hrs = 6         -- rotor yoke height
alpham = 0.88   -- pole ratio
ls = ${model['rlength']*1e3}         -- axial length
mtype = 2       -- magnet shape: 1 = arc, 2 = rectangular

Qm = Q/2        -- slot number in model
Pm = P/2        -- pole number in model

urs = 1000      
urr = 1000      
Br = 1.2        
urm = 1.05       

      
-- coordinate

x = {}
y = {} 
for i=1, 15 do
  x[i]=0
  y[i]=0
end

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

-- nodechains

agnp = 1         -- airgap 
ndt(ag)

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

-- mesh

mesh.con1 = 0.1           

create_mesh_se(Da/2-hj/2,0+hj/2)
create_mesh_se((Da+Di)/4,s/4)
create_mesh_se(Di/2+h1/2,s/4)

-- subregions

def_new_subreg(Da/2-hj/2,0+hj/2,"Stator",11)

-- copy and rotate

mirror_nodechains(x[2],y[2],x[15],y[15])

x0,y0 = pd2c(Di/2-ag/3,0)
x1,y1 = pd2c(Da/2,0)
x2,y2 = pd2c(Da/2,360/Q)
x3,y3 = pd2c(Di/2-ag/3,360/Q)
rotate_copy_nodechains(x0,y0,x1,y1,x2,y2,x3,y3,Qm-1)


----------------------------
------ winding -------------
----------------------------

tauq = 360/Q               
Rq = (Di/2+Da/2-hj)/2      

x11,y11 = pd2c(Rq,tauq-tauq/4)
x12,y12 = pd2c(Rq,tauq+tauq/4)

m.wdg_location   =       1.00    --  Wdg location: -1: intern 1: stator; 2: rotor
m.xcoil_1         =     x11 --   center coordinate of 1. coil side [mm]
m.ycoil_1         =     y11--   center coordinate of 1. coil side [mm]
m.xcoil_2         =      x12 --   center coordinate of 2. coil side [mm]
m.ycoil_2         =      y12--   center coordinate of 2. coil side [mm]





