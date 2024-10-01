
m.magn_rad       = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model['magn_height']*1e3}
m.magn_width      =    ${model['magn_width']*1e3}
m.gap_ma_iron     =    ${model['gap_ma_iron']*1e3}
m.air_triangle    =    ${model.get('air_triangle', 0)}
m.iron_height     =    ${model['iron_height']*1e3}
m.magn_rem        =    ${model['magn_rem']}
m.shaft_rad      =     ${model['condshaft_r']*1e3}
m.magn_ori        =    ${model['magn_ori']}
m.bridge_height   =    ${model['bridge_height']*1e3}
m.bridge_width    =    ${model['bridge_width']*1e3}
m.iron_shape      =    ${model['iron_shape']*1e3}

m.zeroangl        =    ${model.get('zeroangle',0)}

m.mcvkey_yoke     =   mcvkey_yoke
m.nodedist        =   ${model.get('nodedist',1)}

 pre_models("Magnet in Iron")

%if model.get('mcvkey_magnet', ''):
gamma = 0
for i = 0, m.npols_gen-1 do
    alfa = (2*i+1)*180/m.num_poles
    x0, y0 = pd2c((m.magn_rad-m.magn_height/2)*math.cos(math.pi/m.num_poles), alfa)
    if i < 2 and m.npols_gen > 1 then
        delete_sreg(x0, y0)
    end
    if m.orient == "cartaniso" or m.orient == "cartiso" then
        gamma = alfa
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, "red", m.mcvkey_magnet, gamma, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, "green", m.mcvkey_magnet, gamma-180, m.orient, m.magncond, m.rlen)
    end
end
%endif
%if model.get('thcond', 0) and model.get('thcap', 0):
beta = math.pi/m.num_poles
xrb,yrb = pr2c(m.yoke_rad+0.1, beta)  -- rotor lamination
thcond = ${model['thcond']}
thcap = ${model['thcap']}
density = ${model.get('density')*1e3}
def_mat_therm(xrb,yrb,'blue',density,thcond,thcap,1)

%if model.get('thcond_shaft', 0) and model.get('thcap_shaft', 0):
if m.shaft_rad < m.yoke_rad then
   thcond = ${model['thcond_shaft']}
   thcap = ${model['thcap_shaft']}
   density = ${model['spmaweight_shaft']*1e3}
   def_mat_therm(m.yoke_rad/2,0.1,'blue',density,thcond,thcap,1) -- Shaft
end
%endif
rm = m.magn_rad*math.cos(beta)
density = ${model['spmaweight_magnet']*1e3}
thcond = ${model.get('thcond_magnet', 8)}
thcap = ${model.get('thcap_magnet', 440)}
for i = 1,m.npols_gen do -- Magnete
  alfa = (2*i-1) * beta
  xmx,ymx = pr2c(rm,alfa)
  def_mat_therm(xmx,ymx,darkgreen-i%2,density,thcond,thcap,1)
end

if m.air_triangle > 0 then
  rs = math.sqrt(rm^2 + m.magn_width^2/4)
  gam = math.atan(m.magn_width/2, rm)
  for i = 1,m.npols_gen do -- air triangle
    alfa = (2*i-1) * beta
    xsx,ysx = pr2c(rs,alfa + gam + 1e-2)
    def_mat_therm(xsx,ysx,skyblue,1.12,0.026,1007,1)
    xsx,ysx = pr2c(rs,alfa - gam - 1e-2)
    def_mat_therm(xsx,ysx,skyblue,1.12,0.026,1007,1)
  end
end


-- add air layer (inside) for heat transfer
h = 3.8
beta = 360*m.npols_gen/m.num_poles


x0, y0 = pd2c(m.shaft_rad, m.zeroangl)
x1, y1 = pd2c(m.shaft_rad-h, m.zeroangl)
x2, y2 = pd2c(m.shaft_rad-h, beta+m.zeroangl)
x3, y3 = pd2c(m.shaft_rad, beta+m.zeroangl)
nc_line(x0, y0, x1, y1, 0)
nc_circle(x1, y1, x2, y2, 0)
nc_line(x2, y2, x3, y3, 0)
x0, y0 = pd2c(m.shaft_rad-h/2, beta/2+m.zeroangl)
create_mesh_se(x0, y0)

%endif
