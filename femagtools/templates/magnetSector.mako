
m.magn_rad       = da2/2
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2


m.magn_height     =    ${model.get(['magnet','magnetSector', 'magn_height'])*1e3}
m.magn_width      =    ${model.get(['magnet','magnetSector', 'magn_width_pct'])*100}
m.condshaft_r     =    ${model.get(['magnet','magnetSector', 'condshaft_r'])*1e3}
m.magn_num        =    ${model.get(['magnet','magnetSector', 'magn_num'])}
m.magn_perm       =    ${model.get(['magnet','magnetSector', 'magn_rfe'])*1e3}
m.magn_l          =    ${model.get(['magnet','magnetSector', 'magn_len'])*100}
m.magn_ori        =    ${model.get(['magnet','magnetSector', 'magn_ori'])}
m.magn_type       =    ${model.get(['magnet','magnetSector', 'magn_type'])}
m.magn_shape      =    ${model.get(['magnet','magnetSector', 'magn_shape'])*1e3}
m.br_height       =    ${model.get(['magnet','magnetSector', 'bridge_height'])*1e3}
m.br_width        =    ${model.get(['magnet','magnetSector', 'bridge_width'])*1e3}

m.zeroangl        =     0.0
m.cond_shaft      =     0.000
m.mcvkey_yoke     =     mcvkey_yoke
m.mcvkey_mshaft   =     mcvkey_shaft

m.nodedist        =   ${model.magnet.get('nodedist',1)}

 pre_models("Magnet-Sector")

%if model.get_mcvkey_magnet():
gamma = 0
for i = 0, m.npols_gen-1 do
    alfa = (2*i+1)*180/m.num_poles
    x0, y0 = pd2c(m.magn_rad - m.magn_height/2, alfa)
    if i < 2 and m.npols_gen > 1 then
        delete_sreg(x0, y0)
    end
    if m.orient == mcartaniso or m.orient == mcartiso then
        gamma = alfa
    end
    if i % 2 == 0 then
	  def_mat_pm_nlin(x0, y0, red, m.mcvkey_magnet, gamma, m.orient, m.magncond, m.rlen)
    else
        def_mat_pm_nlin(x0, y0, green, m.mcvkey_magnet, gamma-180, m.orient, m.magncond, m.rlen)
    end
end
%endif
