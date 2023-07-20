-- EXPERIMENTAL
-- Caution: this will likely be modified
--
slots_gen = get_mod_data("num_slots")
num_slots = get_dev_data("num_slots")
state_of_problem("therm_static")

beta = 360*slots_gen/num_slots
dy1 = ${model['outer_diam']*1e3}
dy2 = ${model['inner_diam']*1e3}
m.zeroangl = 0
-- heat transfer outside
heat_transfer_coefficient = ${model['heat_transfer_coefficient'][0]}
area_factor = 1
x,y = pd2c(dy1/2+0.05, beta/2+m.zeroangl)
def_heat_transfer(x,y,yellow,heat_transfer_coefficient, area_factor)
-- heat transfer inside
heat_transfer_coefficient = ${model['heat_transfer_coefficient'][1]}
area_factor = 1
x,y = pd2c(dy2/2-0.05, beta/2+m.zeroangl)
def_heat_transfer(x,y,yellow,heat_transfer_coefficient, area_factor)

---------------------------------------------
--  import losses
import_losses_from_femag_dc()
-- overwrite magnet losses (IALH)
%if not isinstance(model.get('magnet_loss_th', 0), list):
printf('magnet losses from B2 Method')
%else:
printf('magnet losses from IAlH Method')
%for i in model['magnet_loss_th']: 
elkeys = get_spel_data("elkeys", ${i[0]})
x, y = get_elem_data("xycp", elkeys[1])
def_losses ( x,y,'red', ${i[1]})
%endfor
%endif

color_gradation_th(0,0,tot,Losses,0,0,model.."_losses.svg")

---------------------------------------------
calc_therm_field()
color_gradation_th(0,0,tot,Temp,0.0,0,model.."_temp.svg")
