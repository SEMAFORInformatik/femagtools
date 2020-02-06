
m.magn_rad       = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model['magn_height']*1e3}
m.magn_width      =    ${model['magn_width']*1e3}
m.gap_ma_iron     =    ${model['gap_ma_iron']*1e3}
m.iron_shape      =    ${model['iron_shape']*1e3}
m.air_space_h     =    ${model['air_space_h']*1e3}
m.iron_bfe        =    ${model['iron_bfe']*1e3}
m.magn_di_ra      =    ${model['magn_di_ra']*1e3}
m.corner_r        =    ${model['corner_r']*1e3}
m.air_sp_ori      =    ${model['air_sp_ori']}
m.magn_ori        =    ${model['magn_ori']}
m.magn_num        =    ${model['magn_num']}
m.cond_shaft      =    0.0 -- ignored
m.zeroangl        =     0.0

m.mcvkey_yoke     =   mcvkey_yoke
m.mcvkey_mshaft    =  mcvkey_shaft
m.nodedist        =   ${model.get('nodedist',1)}

 pre_models("Magnet Iron 4")

%if model.get('mcvkey_magnet', ''):
gamma = 0
for i = 0, m.npols_gen-1 do
    alfa = (2*i+1)*180/m.num_poles
    x0, y0 = pd2c(m.magn_rad - m.magn_height/2 - m.magn_di_ra, alfa)
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
