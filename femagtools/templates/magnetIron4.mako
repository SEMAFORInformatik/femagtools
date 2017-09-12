
m.magn_rad       = da2/2
m.yoke_rad       = dy2/2

m.magn_height     =    ${model.get(['magnet','magnetIron4', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetIron4', 'magn_width'])*1e3}
m.gap_ma_iron     =    ${model.get(['magnet','magnetIron4', 'gap_ma_iron'])*1e3}
m.iron_shape      =    ${model.get(['magnet','magnetIron4', 'iron_shape'])*1e3}
m.air_space_h     =    ${model.get(['magnet','magnetIron4', 'air_space_h'])*1e3}
m.iron_bfe        =    ${model.get(['magnet','magnetIron4', 'iron_bfe'])*1e3}
m.magn_di_ra      =    ${model.get(['magnet','magnetIron4', 'magn_di_ra'])*1e3}
m.corner_r        =    ${model.get(['magnet','magnetIron4', 'corner_r'])*1e3}
m.air_sp_ori      =    ${model.get(['magnet','magnetIron4', 'air_sp_ori'])}
m.magn_ori        =    ${model.get(['magnet','magnetIron4', 'magn_ori'])}
m.magn_num        =    ${model.get(['magnet','magnetIron4', 'magn_num'])}
m.cond_shaft      =    0.0 -- ignored
m.zeroangl        =     0.0

m.mcvkey_yoke     =   mcvkey_yoke
m.mcvkey_mshaft    =  mcvkey_shaft
m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet Iron 4")

%if model.get_mcvkey_magnet():
for i = 0, m.npols_gen-1 do
    x0, y0 = pr2c(m.magn_rad - m.magn_height/2 - m.magn_di_ra,
                  (2*i+1)*math.pi/m.num_poles)
    if i % 2 == 0 then
        def_mat_pm_nlin(x0, y0, red, m.mcvkey_magnet, 0, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, green, m.mcvkey_magnet, 180, m.orient, m.magncond, m.rlen)
    end
end
%endif
