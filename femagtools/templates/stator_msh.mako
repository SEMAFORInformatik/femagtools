

read_msh("${model['name']}")
x1, y1 = pd2c(dy1/2, 360/m.tot_num_slot/2)
x2, y2 = pd2c(da1/2, 360/m.tot_num_slot/2)
mirror_nodechains(x1, y1, x2, y2)

x1, y1 = pd2c(da1/2, 0)
x2, y2 = pd2c(dy1/2, 0)
x3, y3 = pd2c(dy1/2, 360/m.tot_num_slot)
x4, y4 = pd2c(da1/2, 360/m.tot_num_slot)
rotate_copy_nodechains(x1, y1, x2, y2, x3, y3, x4, y4, m.num_sl_gen-1)

-- airgap meshing

x1, y1 = pd2c(da1/2, 0)
x2, y2 = pd2c(da1/2 - ag/3, 0)
x3, y3 = pd2c(da1/2 - ag/3, 360/m.tot_num_slot*m.num_sl_gen)
x4, y4 = pd2c(da1/2, 360/m.tot_num_slot*m.num_sl_gen)

nc_line(x1, y1, x2, y2, 2)
nc_line(x3, y3, x4, y4, 2)
nc_circle_m(x2, y2, x3, y3, 0, 0, num_agnodes+1)
x1, y1 = pd2c(da1/2 - ag/6, 360/m.tot_num_slot*m.num_sl_gen/2)
create_mesh_se(x1, y1)

% for y in model['yoke']:
  x, y = ${model[y][0]*1e3}, ${model[y][1]*1e3}
  if( mcvkey_yoke ~= 'dummy' ) then
    def_mat_fm_nlin(x, y, "blue", mcvkey_yoke, 100) -- ${y}
  else
    def_mat_fm(x, y, "blue", 1000, 100) -- ${y}
  end
% endfor
% for t in model['teeth']:
  x, y = ${model[t][0]*1e3}, ${model[t][1]*1e3}
  if( mcvkey_yoke ~= 'dummy' ) then
    def_mat_fm_nlin(x, y, "blue", mcvkey_teeth, 100) -- ${t}
  else
    def_mat_fm_nlin(x, y, "blue", 1000, 100) -- ${t}
  end
% endfor
% for a in model['air']:
  x, y = ${model[a][0]*1e3}, ${model[a][1]*1e3}
  def_mat_air(x, y)  -- ${a}
% endfor
% if 'wdg' in model:
m.xcoil_1, m.ycoil_1 = ${model['wdg'][0]*1e3}, ${model['wdg'][1]*1e3}
delete_sreg(m.xcoil_1, m.ycoil_1)

r, phi = c2pd(m.xcoil_1, m.ycoil_1)
for n = 1, m.num_sl_gen do
  s = tostring(n)
  x, y = pd2c( r, (n-1)*360/m.tot_num_slot + phi)
  def_new_sreg(x, y, s, 'white')
  x, y = pd2c(r, n*360/m.tot_num_slot - phi)
  add_to_sreg(x, y, s)
end
% else:
m.xcoil_1, m.ycoil_1 = pd2c(da1/2 + (dy1-da1)/8, 360/m.tot_num_slot/2-0.1)
create_mesh()
% endif
