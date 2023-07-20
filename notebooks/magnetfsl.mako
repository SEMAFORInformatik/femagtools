
------------------------------
----- Rotor ------------------
------------------------------
bm = ${model['bm']*1e3}
hm = ${model['hm']*1e3}

-- coordinate

x[1],y[1] = pd2c(Di/2-ag*2/3,0)
x[2],y[2] = pd2c(Di/2-ag*2/3,360/P)
x[3],y[3] = pd2c(Di/2-ag,0)
x[4],y[4] = pd2c(Di/2-ag,360/P)
x[5],y[5] = pd2c(Di/2-ag-hm,0)
x[6],y[6] = pd2c(Di/2-ag-hm,360/P)
x[7],y[7] = pd2c(Di/2-ag-hm-hrs,0)
x[8],y[8] = pd2c(Di/2-ag-hm-hrs,360/P)
x[8],y[8] = pd2c(Di/2-ag-hm-hrs,360/P)

if (mtype==1) then    -- arc
  x[9],y[9]   = pd2c(Di/2-ag,360/P*(1-alpham)*0.5)
  x[10],y[10] = pd2c(Di/2-ag,360/P*(1+alpham)*0.5)
  x[11],y[11] = pd2c(Di/2-ag-hm,360/P*(1-alpham)*0.5)
  x[12],y[12] = pd2c(Di/2-ag-hm,360/P*(1+alpham)*0.5)
else  -- rectangular
  R = Di/2-ag-hm
  alpha = math.asin(bm/(2*R))
  x[11],y[11] = pd2c(R,360/P*0.5-alpha/math.pi*180)
  x[12],y[12] = pd2c(R,360/P*0.5+alpha/math.pi*180)
  R = Di/2-ag
  alpha = math.asin(bm/(2*R))
  x[9],y[9] = pd2c(R,360/P*0.5-alpha/math.pi*180)
  x[10],y[10] = pd2c(R,360/P*0.5+alpha/math.pi*180)
end

-- nodechains

ndt(0.5)

nc_circle(x[1],y[1],x[2],y[2],360/P/agnp+1)
nc_circle(x[3],y[3],x[9],y[9],0)
nc_circle(x[9],y[9],x[10],y[10],0)
nc_circle(x[10],y[10],x[4],y[4],0)
nc_circle(x[5],y[5],x[11],y[11],0)
if (mtype==1) then
  nc_circle(x[11],y[11],x[12],y[12],0)
else
  nc_line(x[11],y[11],x[12],y[12],0)
  nc_line(x[9],y[9],x[10],y[10],0)
end
nc_circle(x[12],y[12],x[6],y[6],0)
nc_circle(x[7],y[7],x[8],y[8],0)

nc_line(x[7],y[7],x[5],y[5],0)
nc_line_cont(x[3],y[3],0)
nc_line_cont(x[1],y[1],0)
nc_line(x[8],y[8],x[6],y[6],0)
nc_line_cont(x[4],y[4],0)
nc_line_cont(x[2],y[2],0)
nc_line(x[11],y[11],x[9],y[9],0)
nc_line(x[12],y[12],x[10],y[10],0)

-- mesh

create_mesh_se(pd2c(Di/2-ag*5/6,180/P))
create_mesh_se(pd2c(Di/2-ag-hm/2,180/P))
create_mesh_se(pd2c(Di/2-ag-hm-hrs/2,180/P))
create_mesh_se(pd2c(Di/2-ag-hm/2,360/P*(1-alpham)*0.25))
create_mesh_se(pd2c(Di/2-ag-hm/2,360/P*(3+alpham)*0.25))

if (mtype==2) then
  R = Di/2-ag-hm
  alpha = math.asin(bm/(2*R))
  Rx = R*math.cos(alpha)
  create_mesh_se(pd2c(Di/2-ag-(R-Rx)/2,360/P*0.5))
end

-- subregions

def_new_subreg(Di/2-ag-hm-hrs/2,ag,"RTYK",11)

-- rotate and copy

rotate_copy_nodechains(x[7],y[7],x[1],y[1],x[2],y[2],x[8],y[8],Pm-1)

-------------------------------
---- Material Properties ----
-------------------------------

-- stator and magnet back iron

def_mat_fm(Da/2-hj/2,ag,urs,100)
def_mat_fm(Di/2-ag-hm-hrs/2,ag,urr,100)

-- magnet

if (mtype==2) then
  for i=0, Pm/2 do
  	alpha = 360/P*(2*i+1)-180/P
    x,y = pd2c(Di/2-ag-hm/2,alpha)
    def_mat_pm(x,y,red,Br,urm,alpha,m.parallel,100)
  end
  for i=1, Pm/2 do
  	alpha = 360/P*2*i-180/P
    x,y = pd2c(Di/2-ag-hm/2,alpha)
    def_mat_pm(x,y,green,Br,urm,alpha+180,m.parallel,100)
  end
else
  for i=0, Pm/2 do
    x,y = pd2c(Di/2-ag-hm/2,360/P*(2*i+1)-180/P)
    def_mat_pm(x,y,red,Br,urm,0,m.radial,100)
  end
  for i=1, Pm/2 do
    x,y = pd2c(Di/2-ag-hm/2,360/P*2*i-180/P)
    def_mat_pm(x,y,green,Br,urm,180,m.radial,100)
  end
end

-- Magnet Data
 
m.remanenc       =         Br     --   Remanence  Br  (Ref:20 Degree C)  [T]   
m.relperm        =        urm     --   Rel. Permeability muer                  
m.spmaweight     =       7.60     --   Specific Weight Magnets       [gr/cm3]  
m.temcoefbr      =      -0.09     --   Temperature Coefficient for Br   [%/K]  
m.temcoefhc      =     -0.100     --   Temperature Coefficient for Hc   [%/K]  
m.magntemp       =       20.0     --   Magnet Temperature          [Degree C]  
m.magncond       =    0.625e6     --   Magnet el. conductivity      [1/Ohm m]  
m.magsegwid      =         bm     --   Magnet segment width              [mm]  
m.magseglen      =         ls     --   Magnet segment length z-direction [mm]  
 
pre_models("Magnet-data")
