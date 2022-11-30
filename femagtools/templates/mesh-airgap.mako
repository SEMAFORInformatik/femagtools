
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
end
