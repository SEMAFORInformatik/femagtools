-- Winding (ROTOR_KS2)

Q1 = m.tot_num_sl
Q2 = ${model['num_slots']}
m.tot_num_sl      = Q2
% if model.get('num_slots_gen', 0):
m.num_sl_gen  =   ${model['num_slots_gen']}
% else:
m.num_sl_gen  =   Q2 * m.npols_gen/m.num_poles
% endif



m.zeroangl        = ${model.get('zeroangl',0)}

m.mcvkey_yoke = mcvkey_yoke


dphi=1e-2

num_wires=1
dir = {'wi', 'wo'}

-- set subregions 'Bar'
r=da2/2 - m.slot_height/2
phi=180/Q2+m.zeroangl
xcoil, ycoil = pd2c(r,phi-dphi)
--delete_sreg(xcoil, ycoil)
--def_new_wdg(xcoil, ycoil, yellow, "Exc", num_wires, 10.0, dir[1])
def_currdens_se(xcoil, ycoil, 5)

for i=2, m.num_sl_gen do
  phi=(2*i-1)*180/Q2 + m.zeroangl
  xcoil, ycoil = pd2c(r,phi)
  --delete_sreg(xcoil, ycoil)
  pole_per_encoche = Q2 / m.num_poles
  pole_num = i // pole_per_encoche

  if pole_num % 2 == 0
	then def_currdens_se(xcoil, ycoil, 5) --add_to_wdg(xcoil, ycoil, wsamekey, dir[1], 'wser')
  else
     def_currdens_se(xcoil, ycoil, -5) --add_to_wdg(xcoil, ycoil, wsamekey, dir[2], 'wser')
end
end
