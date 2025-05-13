
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
    if(m.npols_gen == m.num_poles) then
      alfa = alfa / 2
    end
    x1, y1 = pr2c(r1, alfa)
    n = math.floor(m.fc_radius*alfa/agndst + 1.5)
    nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)
    if(m.npols_gen == m.num_poles) then
       nc_circle_m(x1, y1, r1, 0, 0.0, 0.0, n)
    end

    r2 = m.fc_radius + ag/6
    x2, y2 = pr2c(r2, alfa)
    nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)
    if(m.npols_gen == m.num_poles) then
       nc_circle_m(x2, y2, r2, 0, 0.0, 0.0, n)
    end

    if inner_da_start == nil then
      inner_da_start = da2/2
    end
    x1, y1 = pr2c(inner_da_start, 0.0)
    nc_line(x1, y1, r1, 0.0, 0.0)

    if outer_da_start == nil then
      outer_da_start = da1/2
    end
    x2, y2 = pr2c(outer_da_start, 0.0)
    nc_line(r2, 0.0, x2, y2, 0.0)

    if m.tot_num_slot > m.num_sl_gen then
      if inner_da_end == nil then
        inner_da_end = inner_da_start
      end
      x3, y3 = pr2c(inner_da_end, alfa)
      x4, y4 = pr2c(r1, alfa)
      nc_line(x3, y3, x4, y4, 0, 0)

      if outer_da_end == nil then
        outer_da_end = outer_da_start
      end
      x3, y3 = pr2c(outer_da_end, alfa)
      x4, y4 = pr2c(r2, alfa)
      nc_line(x3, y3, x4, y4, 0, 0)
    end

    x0, y0 = pr2c(r1-ag/6, alfa/2)
    create_mesh_se(x0, y0)
    x0, y0 = pr2c(r2+ag/6, alfa/2)
    create_mesh_se(x0, y0)
% else:
  -- airgap nodechains for axial flux
  x1, y1 = 0, -ag/2 -- airgap center
  x2, y2 = m.num_slots*(m.tooth_width+m.slot_width), y1
  nc_line(x1, -ag/3, x2, -ag/3, num_agnodes+1)
  nc_line(x1, -ag/3, x1, 0, 1)
  nc_line(x2, -ag/3, x2, 0, 1)
  create_mesh_se((x1+x2)/2, -ag/6)

  nc_line(x1, -2*ag/3, x2, -2*ag/3, num_agnodes+1)
  nc_line(x1, -2*ag/3, x1, -ag, 1)
  nc_line(x2, -2*ag/3, x2, -ag, 1)
  create_mesh_se((x1+x2)/2, -5*ag/6)

  nc_line(x1, -ag/3, x1, -2*ag/3, 1)
  nc_line(x2, -ag/3, x2, -2*ag/3, 1)
  create_mesh_se((x1+x2)/2, -ag/2)

  if m.model_type == 'S2R1_all' then
    x1,y1 = 0, -ag -m.magn_height/2
    x2,y2 = m.npols_gen*m.pole_width, y1
    mirror_nodechains(x1,y1, x2,y2)
-- TODO: only needed for geom preview
    x,y = (P22.x+P32.x)/2, -ag - 3*m.magn_height/4
    for i = 1,m.npols_gen do
      if i%2 == 0 then
        def_mat_pm(x,y, 'red', m.remanenc, m.relperm, 90, "parallel", 0, m.rlen)
      else
        def_mat_pm(x,y, 'green', m.remanenc, m.relperm, 270, "parallel", 0, m.rlen)
      end
      x = x + m.pole_width
    end
  elseif m.model_type == 'S1R2_all' then
    x1,y1 = m.npols_gen*m.pole_width, m.slot_height/2
    x2,y2 = 0, y1
    mirror_nodechains(x1,y1, x2,y2)
-- TODO: only needed for geom preview
    x,y = (P22.x+P32.x)/2, m.slot_height - (P21.y+P22.y)/2
    for i = 1,m.npols_gen do
      if i%2 == 0 then
        def_mat_pm(x,y, 'red', m.remanenc, m.relperm, 90, "parallel", 0, m.rlen)
      else
        def_mat_pm(x,y, 'green', m.remanenc, m.relperm, 270, "parallel", 0, m.rlen)
      end
      x = x + m.pole_width
    end
  end

  --  set boundary conditions
  del_bcond()
  if m.st_yoke_height > 0 then
    sh = m.slot_height
  else
    sh = m.slot_height/2
  end
  x1,y1 = 0, sh+m.st_yoke_height
  if m.model_type == 'S2R1_all' or m.model_type == 'S1R2_all' then
    x2,y2 = x1, -ag - m.magn_height -m.yoke_height - y1
  else
    x2,y2 = x1, -ag - m.magn_height -m.yoke_height
  end
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
