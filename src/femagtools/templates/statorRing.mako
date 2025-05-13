-- stator ring model
-- (used for eddy current simulation of 1 rotor bar)
--

zeroangle = m.zeroangl*math.pi/180

post_models("nodedistance", "ndst" )
agndst=ndst[1]*1e3
ndt(agndst)

dnd = ndst[2]/180*math.pi/2
xw1, yw1 = pr2c(da2/2+2*ag/3, zeroangle + math.pi/Q2- dnd)
xw2, yw2 = pr2c(da2/2+2*ag/3, zeroangle + math.pi/Q2 + dnd)
xw3, yw3 = pr2c(da2/2+ag, zeroangle + math.pi/Q2- dnd)
xw4, yw4 = pr2c(da2/2+ag, zeroangle + math.pi/Q2 + dnd)
nc_line(xw1, yw1, xw2, yw2, 2)
nc_line_cont(xw4, yw4, 2)
nc_line_cont(xw3, yw3, 2)
nc_line_cont(xw1, yw1, 2)
  m.xcoil_1, m.ycoil_1 = pr2c((da1-ag/3)/2, zeroangle + math.pi/Q2)
  create_mesh_se(m.xcoil_1, m.ycoil_1)
  stator=def_new_wdg(m.xcoil_1, m.ycoil_1, "green", "1", 1, 0.0, "wi")

x1, y1 = pr2c(da2/2, zeroangle)
x2, y2 = pr2c(da2/2+ag/3, zeroangle)
x3, y3 = pr2c(da2/2+2*ag/3, zeroangle)
x4, y4 = pr2c(da2/2+ag, zeroangle)
x5, y5 = pr2c(da1/2 + (dy1-da1)/12, zeroangle)
x6, y6 = pr2c(da1/2 + (dy1-da1)/12, zeroangle + 2*math.pi/Q2)
x7, y7 = pr2c(da2/2+ag, zeroangle + 2*math.pi/Q2)
x8, y8 = pr2c(da2/2+2*ag/3, zeroangle + 2*math.pi/Q2)
x9, y9 = pr2c(da2/2+ag/3, zeroangle + 2*math.pi/Q2)
x10, y10 = pr2c(da2/2, zeroangle + 2*math.pi/Q2)

nc_line(x1, y1, x2, y2, 0)
nc_line_cont(x3, y3, 0)
nc_line_cont(x4, y4, 0)
ndt(2.5*agndst)
nc_line_cont(x5, y5, 0)
ndt(4*agndst)
nc_circle_m(x5, y5, x6, y6, 0,0,0)
ndt(2.5*agndst)
nc_line(x6, y6, x7, y7, 0)
ndt(agndst)
nc_line(x7, y7, x8, y8, 0)
nc_line_cont(x9, y9, 0)
nc_line_cont(x10, y10, 0)

nc_circle_m(x2, y2, x9, y9, 0,0,0)
nc_circle_m(x3, y3, xw1, yw1, 0,0,0)
nc_circle_m(xw2, yw2, x8, y8, 0,0,0)
nc_circle_m(x4, y4, xw3, yw3, 0,0,0)
nc_circle_m(xw4, yw4, x7, y7, 0,0,0)

x, y = (x4+x7)/2, (y4+y5)/2
create_mesh_se(x, y)

def_new_sreg(x, y, 'StZa', 'skyblue')
  if mcvkey_yoke ~= 'dummy' then
    def_mat_fm_nlin(x, y, 'blue', mcvkey_yoke, 100)
  else
    def_mat_fm(x, y, ur, 100)
  end
create_mesh()

x0, y0 = pr2c( dy2/2, math.pi*(1/2 - 1/Q2))
x9, y9 = pr2c( dy2/2, math.pi*(1/2 + 1/Q2))
def_bcond_vpo(x0,y0, x4,y4)
def_bcond_vpo(x4,y4, x5,y5)
def_bcond_vpo(x5,y5, x9,y9)
