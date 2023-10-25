
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

m.zeroangl        =    ${model.get('zeroangle',0)}
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
%if model.get('thcond', 0) and model.get('thcap', 0):
beta = math.pi/m.num_poles + m.zeroangl/180*math.pi
xrb,yrb = pr2c(m.yoke_rad+0.1, beta)  -- rotor lamination
thcond = 24 -- ${model['thcond']}
thcap = 480 -- ${model['thcap']}
def_mat_therm(xrb,yrb,'blue',7700,thcond,thcap,1)
--def_mat_therm(m.yoke_rad/2,0.1,'blue',7700,thcond,thcap,1) -- Shaft

rm = m.rotor_rad - m.magn_height/2
thcond = 8
thcap = 440
for i = 1,m.npols_gen do -- Magnets
  alfa = (2*i-1) * beta
  xmx,ymx = pr2c(rm,alfa-beta/2)
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
--[[
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
--]]
%endif