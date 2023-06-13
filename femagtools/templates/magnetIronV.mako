
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model['magn_height']*1e3}
m.magn_width      =    ${model['magn_width']*1e3}
m.magn_angle      =    ${model['magn_angle']}
m.magn_num        =    ${model['magn_num']}
m.iron_hs         =    ${model['iron_hs']*1e3}
m.iron_height     =    ${model['iron_height']*1e3}
m.iron_shape      =    ${model['iron_shape']*1e3}
m.air_triangle    =    ${model.get('air_triangle', 0)}
m.gap_ma_iron     =    ${model['gap_ma_iron']*1e3}
m.magn_rem        =    ${model.get('magn_rem', 'm.remanenc')}
m.shaft_rad       =    ${model['condshaft_r']*1e3}

m.zeroangl        =    ${model.get('zeroangle',0)}

m.mcvkey_yoke     =   mcvkey_yoke

m.nodedist        =   ${model.get('nodedist',1)}

 pre_models("Magnet Iron V")

%if model.get('mcvkey_magnet', ''):
alpha = math.pi/m.num_poles/2

beta = {}
beta[1] = m.magn_angle/360*math.pi - 2*alpha
beta[2] = m.magn_angle/360*math.pi + 2*alpha

r = m.rotor_rad
d = r* (math.sin(4*alpha) - math.tan(beta[2])*math.cos(4*alpha))/
      (math.sin(3*alpha) - math.tan(beta[2])*math.cos(3*alpha)) + m.magn_height

x0, y0 = pr2c(d, alpha)
delete_sreg(x0, y0)
if m.npols_gen > 1 then
   x0, y0 = pr2c(d, 5*alpha)
   delete_sreg(x0, y0)
end

delta = 0
for i = 0, m.npols_gen-1 do
    gamma = (2*i+1)*180/m.num_poles
    x0, y0 = pd2c(d, gamma-alpha/math.pi*180)
    if m.orient == "cartaniso" or m.orient == "cartiso" then
        delta = gamma-m.magn_angle/2+90
    else
        delta = m.magn_angle/2 - 90/m.num_poles
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, "red", m.mcvkey_magnet, delta, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, "green", m.mcvkey_magnet, delta-180, m.orient, m.magncond, m.rlen)
    end
    x0, y0 = pd2c(d, gamma+alpha/math.pi*180)
    if m.orient == "cartaniso" or m.orient == "cartiso" then
       delta = gamma+m.magn_angle/2-90
    else
        delta = -m.magn_angle/2 + 90/m.num_poles
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, "red", m.mcvkey_magnet, delta, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, "green", m.mcvkey_magnet, delta+180, m.orient, m.magncond, m.rlen)
    end
end
%endif

%if model.get('thcond', 0) and model.get('thcap', 0):
beta = math.pi/m.num_poles + m.zeroangl/180*math.pi
xrb,yrb = pr2c(m.yoke_rad+0.1, beta)  -- rotor lamination
thcond = ${model['thcond']}
thcap = ${model['thcap']}
def_mat_therm(xrb,yrb,'blue',7700,thcond,thcap,1)
--def_mat_therm(m.yoke_rad/2,0.1,'blue',7700,thcond,thcap,1) -- Shaft

gam = m.magn_angle/2/180*math.pi
y = m.rotor_rad*math.sin(beta)
x = m.rotor_rad*math.cos(beta)
rr = y/math.sin(gam)
xx = x - math.sqrt(rr^2 - y^2)
rm = xx/(math.cos(beta/2) - math.sin(beta/2)/math.tan(gam))
thcon = 8
thcap = 440
for i = 1,m.npols_gen do -- Magnets
  alfa = (2*i-1) * beta
  xmx,ymx = pr2c(rm,alfa-beta/2)
  def_mat_therm(xmx,ymx,darkgreen-i%2,7500,thcond,thcap,1)
  xmx,ymx = pr2c(rm,alfa+beta/2)
  def_mat_therm(xmx,ymx,darkgreen-i%2,7500,thcond,thcap,1)
end

thcond = 0.026   -- air
thcap = 1007
get_spel_keys("sekeys")      -- Get all subregions of the model
for i=1, #sekeys do
  srkey =  get_spel_data("srkey", sekeys[i])
  if srkey > 0 then
    srname = get_sreg_data("name",srkey)
    if srname == '    ' then
      srkey = 0
    end
  end
  if srkey == 0 then
    elkeys = get_spel_data("elkeys", sekeys[i])
    Ex, Ey = get_elem_data("xycp", elkeys[1])
    r, phi = c2pr(Ex, Ey)
    if r < m.rotor_rad then -- rotor only
      def_mat_therm(Ex,Ey,'skyblue',1.12,thcond,thcap,1)
    end
  end
end

-- add air layer (inside) for heat transfer
h = 3.8
beta = 360*m.npols_gen/m.num_poles
x0, y0 = pd2c(dy2/2, m.zeroangl)
x1, y1 = pd2c(dy2/2-h, m.zeroangl)
x2, y2 = pd2c(dy2/2-h, beta+m.zeroangl)
x3, y3 = pd2c(dy2/2, beta+m.zeroangl)
nc_line(x0, y0, x1, y1, 0)
nc_circle(x1, y1, x2, y2, 0)
nc_line(x2, y2, x3, y3, 0)
x0, y0 = pd2c(dy2/2-h/2, beta/2+m.zeroangl)
create_mesh_se(x0, y0)

%endif