
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
% if model.get(['magnet','magnetIronV', 'magn_rem']):
m.magn_rem        =    ${model.get(['magnet','magnetIronV', 'magn_rem'])}
% else:
m.magn_rem        =    ${model.get('remanenc', 1.2)}
%endif
m.shaft_rad       =     ${model.get(['magnet','magnetIronV', 'condshaft_r'])*1e3}

m.zeroangl        =     0.0

m.mcvkey_yoke     =   mcvkey_yoke

m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet Iron V")

%if model.get_mcvkey_magnet():
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
    if m.orient == mcartaniso or m.orient == mcartiso then
        delta = gamma-m.magn_angle/2+90
    else
        delta = m.magn_angle/2 - 90/m.num_poles
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, red, m.mcvkey_magnet, delta, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, green, m.mcvkey_magnet, delta-180, m.orient, m.magncond, m.rlen)
    end
    x0, y0 = pd2c(d, gamma+alpha/math.pi*180)
    if m.orient == mcartaniso or m.orient == mcartiso then
       delta = gamma+m.magn_angle/2-90
    else
        delta = -m.magn_angle/2 + 90/m.num_poles
    end
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, red, m.mcvkey_magnet, delta, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, green, m.mcvkey_magnet, delta+180, m.orient, m.magncond, m.rlen)
    end
end
%endif
