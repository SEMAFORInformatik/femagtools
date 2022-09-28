-- Winding (ROT_HSM only!)
m.num_wires =  ${model['num_wires']}

a = m.rotor_diam/2 - ag - m.pole_height - m.core_height/2
b = ((m.core_width1 + m.core_width2)/4 + m.pole_width/2)/2
r = math.sqrt(a^2 + b^2)
alpha = math.atan2(b, a)

phi = math.pi/m.num_poles
dir = {'wi', 'wo'}
xcoil, ycoil = pr2c(r, phi - alpha)
def_new_wdg(xcoil, ycoil, yellow, "Exc", m.num_wires, 10.0, dir[1])
xcoil, ycoil = pr2c(r, phi + alpha)
add_to_wdg(xcoil, ycoil, wsamekey, dir[2], 'wser')
for i = 2, m.npols_gen do
  n = (i+1) % 2 + 1
  phi = phi + 2*math.pi/m.num_poles
  xcoil, ycoil = pr2c(r, phi - alpha)
  add_to_wdg(xcoil, ycoil, wsamekey, dir[n], 'wser')
  xcoil, ycoil = pr2c(r, phi + alpha)
  n = i % 2 + 1
  add_to_wdg(xcoil, ycoil, wsamekey, dir[n], 'wser')
end
