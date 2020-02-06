
m.magn_rad       = da2/2
m.rotor_rad      = da2/2
m.yoke_rad       = dy2/2


m.magn_height     =    ${model['magn_height']*1e3}
% if model.get('magn_width_pct', 0):
m.magn_width      =    ${model['magn_width_pct']*100}
% else:
m.magn_width      =    -${model['magn_width']*1e3}
% endif
m.condshaft_r     =    ${model['condshaft_r']*1e3}
m.magn_num        =    ${model['magn_num']}
m.magn_perm       =    ${model['magn_rfe']*1e3}
m.magn_l          =    ${model['magn_len']*100}
m.magn_ori        =    ${model['magn_ori']}
m.magn_type       =    ${model['magn_type']}
m.magn_shape      =    ${model['magn_shape']*1e3}
m.br_height       =    ${model['bridge_height']*1e3}
m.br_width        =    ${model['bridge_width']*1e3}

m.zeroangl        =     0.0
m.cond_shaft      =     0.000
m.mcvkey_yoke     =     mcvkey_yoke
m.mcvkey_mshaft   =     mcvkey_shaft

m.nodedist        =   ${model.get('nodedist',1)}

 pre_models("Magnet-Sector")

%if model.get('magn_ori', 0) == 8:
if m.magncond == nil then
  m.magncond = 6.25e5
end
if m.rlen == nil then
  m.rlen = 100
end

for i = 0, m.npols_gen-1 do
    alfa = (2*i+1)*180/m.num_poles
    x0, y0 = pd2c(m.magn_rad - m.magn_height/2, alfa)
    if i % 2 == 0 then
      color = "red"
      phi = alfa+180/m.num_poles
    else
      phi = -90
      color = "green"
    end
    def_mat_pm(x0, y0, color, m.remanenc, m.relperm,
	               phi, m.radial, m.magncond, m.rlen)
end
%elif model.get('mcvkey_magnet', ''):
gamma = 0
for i = 0, m.npols_gen-1 do
    alfa = (2*i+1)*180/m.num_poles
    x0, y0 = pd2c(m.magn_rad - m.magn_height/2, alfa)
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
