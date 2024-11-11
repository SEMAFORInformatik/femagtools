
m.magn_rad       = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model['magn_height']*1e3}
m.magn_width      =    ${model['magn_width']*1e3}
m.gap_ma_iron     =    ${model['gap_ma_iron']*1e3}
m.iron_bfe        =    ${model['iron_bfe']*1e3}
m.air_space_h     =    ${model['air_space_h']*1e3}
m.corner_r        =    ${model['corner_r']*1e3}
m.air_sp_ori      =    ${model['air_sp_ori']}
m.magn_num        =    ${model['magn_num']}
m.iron_shape      =    ${model['iron_shape']*1e3}
m.air_space_b     =    ${model['air_space_b']*1e3}
m.magn_di_ra      =    ${model['magn_di_ra']*1e3}

m.zeroangl        =   ${model.get('zeroangle',0)}

m.mcvkey_yoke     =   mcvkey_yoke
m.mcvkey_mshaft   =   mcvkey_shaft
m.nodedist        =   ${model.get('nodedist',1)}

 pre_models("Magnet Iron 5")

%if model.get('mcvkey_magnet', ''):
gamma = 0
for i = 0, m.npols_gen-1 do
    alfa = (2*i+1)*180/m.num_poles
    x0, y0 = pd2c(m.magn_rad - m.magn_height/2 - m.magn_di_ra, alfa)
    if i < 2 and m.npols_gen > 1 then
        delete_sreg(x0, y0)
    end
    if m.orient == mcartaniso or m.orient == mcartiso then
        gamma = alfa
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, "red", m.mcvkey_magnet, gamma, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, "green", m.mcvkey_magnet, gamma-180, m.orient, m.magncond, m.rlen)
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