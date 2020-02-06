
m.magn_rad       = da2/2
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model['magn_height']*1e3}
m.iron_bfe        =    ${model['iron_bfe']*1e3}
m.gap_ma_iron     =    ${model['gap_ma_iron']*1e3}
m.air_triangle    =    ${model.get('air_triangle', 0)}
m.iron_height     =    ${model['iron_height']*1e3}
m.gap_ma_rigth    =    ${model['gap_ma_right']*1e3}
m.gap_ma_left     =    ${model['gap_ma_left']*1e3}
m.shaft_rad      =     ${model['condshaft_r']*1e3}
m.magn_num        =    ${model['magn_num']}
m.magn_ori        =    ${model['magn_ori']}
m.iron_shape      =    ${model['iron_shape']*1e3}

m.zeroangl        =     0.0

m.mcvkey_yoke     =   mcvkey_yoke
m.nodedist        =   ${model.get('nodedist',1)}

 pre_models("Magnet Iron 3")
 
%if model.get('mcvkey_magnet', ''):
gamma = math.pi/m.num_poles
alfa = 0
if m.magn_num == 1 then
for i = 0, m.npols_gen-1 do
    x0, y0 = pr2c((m.magn_rad - m.iron_shape)*math.cos(gamma),
                  (2*i+1)*gamma)
    if i < 2 and m.npols_gen > 1 then
        delete_sreg(x0, y0)
    end
    if m.orient == mcartaniso or m.orient == mcartiso then
        alfa = (2*i+1)*gamma*180/math.pi
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, "red", m.mcvkey_magnet, alfa, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, "green", m.mcvkey_magnet, alfa-180, m.orient, m.magncond, m.rlen)
    end
end
else
r = (m.magn_rad - m.iron_shape)*math.cos(gamma/2) - m.magn_height/2
for i = 0, m.npols_gen-1 do
    alfa = (2*i+1)*gamma
    x0, y0 = pr2c(r, (alfa + gamma/2)
    x1, y1 = pr2c(r, (alfa - gamma/2)
    if i < 2 then
        delete_sreg(x0, y0)
	delete_sreg(x1, y1)
    end
    if m.orient == "polaniso" or m.orient == "poliso" then
        alfa = 0
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, "red", m.mcvkey_magnet, (alfa+gamma/2)*180/math.pi, m.orient, m.magncond, m.rlen)
        def_mat_pm_nlin(x1, y1, "red", m.mcvkey_magnet, (alfa-gamma/2)*180/math.pi, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, "green", m.mcvkey_magnet, (alfa+gamma/2)*180/math.pi-180, m.orient, m.magncond, m.rlen)
        def_mat_pm_nlin(x1, y1, "green", m.mcvkey_magnet, (alfa-gamma/2)*180/math.pi-180, m.orient, m.magncond, m.rlen)
    end
end
end

%endif
