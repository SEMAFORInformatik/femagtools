--[[
  Ring rotor (either laminated or pure air)
--]]

x1, y1 = pd2c(da2/2, m.zeroangl)
x2, y2 = pd2c(dy2/2, m.zeroangl)
x3, y3 = pd2c(dy2/2, m.zeroangl + 360/m.num_poles)
x4, y4 = pd2c(da2/2, m.zeroangl + 360/m.num_poles)
ndt(2*agndst)
nc_circle(x1, y1, x4, y4, 0)
ndt(8*agndst)
nc_line(x1, y1, x2, y2, 0)
nc_circle(x2, y2, x3, y3, 0)
nc_line(x3, y3, x4, y4,0)
xm, ym = pd2c((dy2+da2)/4, 90/m.num_poles)
create_mesh_se(xm, ym)
if mcvkey_yoke ~= 'dummy' then
  def_mat_fm_nlin(xm, ym, 'blue', m.mcvkey_yoke, 100)
else
  def_mat_air(xm, ym)
end
if(m.npols_gen>1) then
  rotate_copy_nodechains(x2, y2, x1, y1, x4, y4, x3, y3, m.npols_gen-1)
end
