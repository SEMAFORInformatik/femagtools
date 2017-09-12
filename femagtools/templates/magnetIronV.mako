
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model.get(['magnet','magnetIronV', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetIronV', 'magn_width'])*1e3}
m.magn_angle      =    ${model.get(['magnet','magnetIronV', 'magn_angle'])}
m.magn_num        =    ${model.get(['magnet','magnetIronV', 'magn_num'])}
m.iron_hs         =    ${model.get(['magnet','magnetIronV', 'iron_hs'])*1e3}
m.iron_height     =    ${model.get(['magnet','magnetIronV', 'iron_height'])*1e3}
m.iron_shape      =    ${model.get(['magnet','magnetIronV', 'iron_shape'])*1e3}
m.air_triangle    =    ${model.get(['magnet','magnetIronV', 'air_triangle'])}
m.gap_ma_iron     =    ${model.get(['magnet','magnetIronV', 'gap_ma_iron'])*1e3}
m.magn_rem         =    ${model.get(['magnet','magnetIronV', 'magn_rem'])}
m.shaft_rad       =     ${model.get(['magnet','magnetIronV', 'condshaft_r'])*1e3}

m.zeroangl        =     0.0

m.mcvkey_yoke     =   mcvkey_yoke

m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet Iron V")

%if model.get_mcvkey_magnet():
alpha = {}
alpha[1] = math.pi/m.num_poles/2
alpha[2] = 3*alpha[1]
beta = {}
beta[1] = m.magn_angle/360*math.pi - 2*alpha[1]
beta[2] = m.magn_angle/360*math.pi + 2*alpha[1]

r = m.rotor_rad 
d = r* (math.sin(4*alpha[1]) - math.tan(beta[2])*math.cos(4*alpha[1]))/
      (math.sin(3*alpha[1]) - math.tan(beta[2])*math.cos(3*alpha[1])) + m.magn_height

for i = 0, m.npols_gen-1 do
    gamma = 2*i*math.pi/m.num_poles
    x0, y0 = pr2c(d, alpha[1]+gamma)
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, red, m.mcvkey_magnet, 0, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, green, m.mcvkey_magnet, 180, m.orient, m.magncond, m.rlen)
    end
    x0, y0 = pr2c(d, alpha[2]+gamma)
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, red, m.mcvkey_magnet, 0, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, green, m.mcvkey_magnet, 180, m.orient, m.magncond, m.rlen)
    end
end
%endif
