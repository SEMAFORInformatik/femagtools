hm = ${model['magn_height']}*1e3    -- height of magnets [mm]
hrs = (da2-dy2)/2 - hm

P=m.num_poles
Pm=P*m.num_sl_gen/m.tot_num_sl

r=da2/2

agnp = m.nodedist


x = {}
y = {} 

-- Berechnung der Koordinaten

x[1],y[1] = pd2c(r,0)
x[2],y[2] = pd2c(r,360/P)
x[3],y[3] = pd2c(r-hm,0)
x[4],y[4] = pd2c(r-hm,360/P)
x[5],y[5] = pd2c(r-hm-hrs,0)
x[6],y[6] = pd2c(r-hm-hrs,360/P)

-- node chains

nc_circle(x[1],y[1],x[2],y[2],360/P/agnp+1)
nc_circle(x[3],y[3],x[4],y[4],0)
nc_circle(x[5],y[5],x[6],y[6],0)

nc_line(x[5],y[5],x[3],y[3],0)
nc_line_cont(x[1],y[1],0)
nc_line(x[6],y[6],x[4],y[4],0)
nc_line_cont(x[2],y[2],0)

nc_line(x[1], y[1], x[1]+ag/3, y[1], 0)
x[7], y[7] = pd2c(r+ag/3, 360/P)
nc_circle(x[1]+ag/3,y[1],x[7], y[7],360/P/agnp+1)
nc_line(x[2], y[2], x[7], y[7], 0)

-- meshing

create_mesh_se(r-hm/2,1.0)
create_mesh_se(r-hm-hrs/2,1.0)
x0, y0 = pd2c(r+ag/6, 180/P)
create_mesh_se(x0, y0)

-- Definition of subregions

def_new_subreg(r-hm-hrs/2,1.0,"RotorYoke",11)

-- rotate

rotate_copy_nodechains(x[5],y[5],x[1]+ag/3,y[1],
                       x[7],y[7],x[6],y[6],Pm-1)

-------------------------------
---- Material properties ----
-------------------------------

--  PM yoke
if mcvkey_yoke ~= 'dummy' then
  m.rlength = 100
  def_mat_fm_nlin(r-hm-hrs/2, 1.0, blue, mcvkey_yoke, m.rlength)
else
  urr = 1000
  def_mat_fm(r-hm-hrs/2,1.0,urr,100)
end

-- Permanent magnets
for i=0, Pm-1 do
  x,y = pd2c(da2/2-hm/2,360/P*i+180/P)
  if ( i % 2 == 0 ) then
    def_mat_pm(x,y,green, m.remanenc, m.relperm,180,m.radial,100)
  else
    def_mat_pm(x,y,red,m.remanenc, m.relperm,0,m.radial,100)
  end
end


