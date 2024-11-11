
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

%if model.get('thcond', 0) and model.get('thcap', 0) and model.get('density', 0):
if m.shaft_rad == nil then 
    m.shaft_rad = dy2/2
end 
if m.shaft_rad > dy2/2 then 
    m.shaft_rad = dy2/2
end 
beta = math.pi/m.num_poles
rotor_thcond = ${model['thcond']}
rotor_thcap = ${model['thcap']}
rotor_density = ${model.get('density')}

%if model.get('thcond_shaft', 0) and model.get('thcap_shaft', 0):
if m.shaft_rad < m.yoke_rad then
   shaft_thcond = ${model['thcond_shaft']}
   shaft_thcap = ${model['thcap_shaft']}
   shaft_density = ${model['spmaweight_shaft']*1e3}
   r_shaft = (m.shaft_rad + m.yoke_rad)/2
   x0_shaft, y0_shaft = pd2c(r_shaft, beta/2)
end
%endif

if x0_shaft == nil then
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
end

%endif