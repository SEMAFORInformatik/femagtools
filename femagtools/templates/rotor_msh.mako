
read_msh("${model['name']}")

x1, y1 = pd2c(da2/2, 360/m.num_poles/2)
x2, y2 = pd2c(dy2/2, 360/m.num_poles/2)
mirror_nodechains(x1, y1, x2, y2)

-- airgap
x1, y1 = pd2c(da2/2, 0)
x2, y2 = pd2c(da2/2 + ag/3, 0)
x3, y3 = pd2c(da2/2 + ag/3, 360/m.tot_num_slot*m.num_sl_gen)
x4, y4 = pd2c(da2/2, 360/m.tot_num_slot*m.num_sl_gen)

nc_line(x1, y1, x2, y2, 2)
nc_line(x3, y3, x4, y4, 2)
nc_circle_m(x2, y2, x3, y3, 0, 0, num_agnodes+1)
x1, y1 = pd2c(da2/2 + ag/6, 360/m.tot_num_slot*m.num_sl_gen/2)
create_mesh_se(x1, y1)

% for y in model['yoke']:
  x, y = ${model[y][0]*1e3}, ${model[y][1]*1e3}
  def_mat_fm_nlin(x, y, "blue", mcvkey_yoke, 100) -- ${y}
% endfor
% for a in model['air']:
  x, y = ${model[a][0]*1e3}, ${model[a][1]*1e3}
  def_mat_air(x, y)  -- ${a}
% endfor

-- Magnet
 color = 'red'
 if( m.magncond == nil ) then
   m.magncond = 625000.0
 end
 m.rlen = 100
% for i, m in enumerate(model['mag']['sreg']):
  x, y = ${model[m][0]*1e3}, ${model[m][1]*1e3}
  delete_sreg(x, y)
  alfa = ${model['mag']['axis'][i]-90}  -- ${m}
  def_mat_pm(x, y, color, m.remanenc, m.relperm,
	               alfa, m.parallel, m.magncond, m.rlen)
  r, phi = c2pd(x, y)
  x, y = pd2c(r, 360/m.num_poles-phi)
  alfa = 360//m.num_poles/2 - alfa
  def_mat_pm(x, y, color, m.remanenc, m.relperm,
	               alfa, m.parallel, m.magncond, m.rlen)

% endfor

