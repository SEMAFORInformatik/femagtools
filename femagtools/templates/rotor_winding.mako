-- Winding (ROT_HSM only!)
m.num_wires =  ${model['num_wires']}

xm = {}
ym = {}
a = m.rotor_diam/2 - ag - m.pole_height - m.core_height/2
b = ((m.core_width1 + m.core_width2)/4 + m.pole_width/2)/2
r = math.sqrt(a^2 + b^2)
alpha = math.atan2(b, a)

xm[1], ym[1] = pr2c(r, math.pi/m.num_poles - alpha)
xm[2], ym[2] = pr2c(r, math.pi/m.num_poles + alpha)
def_new_wdg(xm[1], ym[1], yellow, "Exc", m.num_wires, 1000.0, wi)
add_to_wdg(xm[2], ym[2], wsamekey, wo, wser) 
