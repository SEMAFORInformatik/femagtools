
-- create airgap nodechains and mesh if not already created...

if not airgap_created then
    if num_agnodes ~= nil then
      alfa = 2*math.pi*m.npols_gen/m.num_poles
    end
    if agndst == nil then
      agndst = m.nodedist
    end
    ndt(agndst)
    if alfa == nil then
      alfa = m.npols_gen*2*math.pi/m.num_poles
    end
% if hasattr(model, 'bore_diam'):
    r1 = m.fc_radius - ag/6
    x1, y1 = pr2c(r1, alfa)
    n = math.floor(m.fc_radius*alfa/agndst + 1.5)
    nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)

    r2 = m.fc_radius + ag/6
    x2, y2 = pr2c(r2, alfa)
    nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)

    if inner_da_start == nil then
      inner_da_start = da2/2
    end
    nc_line(r1, 0.0, r2, 0.0, 2)

    if outer_da_start == nil then
      outer_da_start = da1/2
    end
    nc_line(r2, 0.0, outer_da_start, 0.0, 0)

    if m.tot_num_slot > m.num_sl_gen then
      x4, y4 = pr2c(outer_da_start, alfa)
      nc_line(x1, y1, x2, y2, 2)
      nc_line(x4, y4, x2, y2, 0)
    end

    x0, y0 = pr2c(r2-ag/6, alfa/2)
    create_mesh_se(x0, y0)
    x0, y0 = pr2c(r2+ag/6, alfa/2)
    create_mesh_se(x0, y0)
% else:
  -- airgap nodechains for axial flux
  x1, y1 = 0, -ag/2 -- airgap center
  x2, y2 = m.num_slots*(m.tooth_width+m.slot_width), y1
  nc_line(x1, -ag/3, x2, -ag/3, num_agnodes)
  nc_line(x1, -ag/3, x1, 0, 0)
  nc_line(x2, -ag/3, x2, 0, 0)
  create_mesh_se((x1+x2)/2, -ag/6)

  nc_line(x1, -2*ag/3, x2, -2*ag/3, num_agnodes)
  nc_line(x1, -2*ag/3, x1, -ag, 0)
  nc_line(x2, -2*ag/3, x2, -ag, 0)
  create_mesh_se((x1+x2)/2, -5*ag/6)

  nc_line(x1, -ag/3, x1, -2*ag/3, 0)
  nc_line(x2, -ag/3, x2, -2*ag/3, 0)
  create_mesh_se((x1+x2)/2, -ag/2)

  --  set boundary conditions
  del_bcond()
  x1,y1 = 0, m.tooth_width/2+m.slot_width/2
  x2,y2 = x1, -ag - m.magn_height -m.yoke_height
  x3,y3 = m.npols_gen*m.pole_width, y2
  x4,y4 = x3, y1
  if m.npols_gen%2 == 1 then
    def_bcond_only(x1,y1, x2,y2, x3,y3, x4,y4, 3)
  else
    def_bcond_only(x1,y1, x2,y2, x3,y3, x4,y4, 4)
  end
  if (m.model_type ~= "S1R2") then
   def_bcond_vpo(x4,y4, x1,y1)
  end
  if (m.model_type ~= "S2R1") then
   def_bcond_vpo(x2,y2, x3,y3)
  end

% endif
end
