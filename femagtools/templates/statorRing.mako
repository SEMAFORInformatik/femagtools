-- stator ring model
-- (used for eddy current simulation of 1 rotor bar)
--
Q2 = ${int(model['num_slots'])}

post_models("nodedistance", "ndst" )
agndst=ndst[1]*1e3
ndt(agndst)
x1, y1 = pr2c(da2/2, math.pi*(1/2 - 1/Q2))
x2, y2 = pr2c(da2/2+ag/2, math.pi*(1/2 - 1/Q2))
x3, y3 = pr2c(da2/2+ag, math.pi*(1/2 - 1/Q2))
x4, y4 = pr2c(da1/2 + (dy1-da1)/12, math.pi*(1/2 - 1/Q2))
x5, y5 = pr2c(da1/2 + (dy1-da1)/12, math.pi*(1/2 + 1/Q2))
x6, y6 = pr2c(da2/2+ag, math.pi*(1/2 + 1/Q2))
x7, y7 = pr2c(da2/2+ag/2, math.pi*(1/2 + 1/Q2))
x8, y8 = pr2c(da2/2, math.pi*(1/2 + 1/Q2))

nc_line(x1, y1, x2, y2, 0)
nc_line_cont(x3, y3, 0)
nc_line_cont(x4, y4, 0)
nc_circle_m(x4, y4, x5, y5, 0,0,0)
nc_line(x5, y5, x6, y6, 0)
nc_line_cont(x7, y7, 0)
nc_line_cont(x8, y8, 0)
nc_circle_m(x3, y3, x6, y6, 0,0,0)
nc_circle_m(x2, y2, x7, y7, 0,0,0)

x, y = (x3+x6)/2, (y3+y4)/2
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
