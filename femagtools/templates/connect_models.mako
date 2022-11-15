% if model.stator['num_slots'] > model.stator['num_slots_gen']:
connect_models()
% else:
x0, y0 = pr2c(da1/2-ag/6, math.pi/2)
create_mesh_se(x0, y0)
x0, y0 = pr2c(da2/2+ag/6, math.pi/2)
create_mesh_se(x0, y0)

x0, y0 = pr2c(m.fc_radius, math.pi/2)
create_mesh_se(x0, y0)

x0,y0 = pd2c(dy1/2,0)
def_bcond_vpo(x0, 0,-x0,0,0)
def_bcond_vpo(-x0,0,x0,0,0)
x0,y0 = pd2c(dy2/2,0)
def_bcond_vpo(x0, 0,-x0,0,0)
def_bcond_vpo(-x0,0,x0,0,0)
% endif