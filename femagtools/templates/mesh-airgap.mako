

-- airgap
ndt(agndst)
r1 = da2/2 + ag/3
x1, y1 = pr2c(r1, alfa)
n = r1*alfa/agndst + 1
nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)

r2 = da2/2 + 2*ag/3
x2, y2 = pr2c(r2, alfa)
nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)

if inner_max_corner_x == nil then
  inner_max_corner_x = da2/2
end
x1, y1 = inner_max_corner_x, 0.0
nc_line(x1, y1, r1, 0.0, 0.0)

if outer_min_corner_x == nil then
  outer_min_corner_x = da1/2
end
x2, y2 = outer_min_corner_x, 0.0
nc_line(r2, 0.0, x2, y2, 0.0)

x3, y3 = pr2c(x1, alfa)
x4, y4 = pr2c(r1, alfa)
nc_line(x3, y3, x4, y4, 0, 0)

x3, y3 = pr2c(x2, alfa)
x4, y4 = pr2c(r2, alfa)
nc_line(x3, y3, x4, y4, 0, 0)

x0, y0 = pr2c(r1-ag/6, alfa/2)
create_mesh_se(x0, y0)
x0, y0 = pr2c(r2+ag/6, alfa/2)
create_mesh_se(x0, y0)
