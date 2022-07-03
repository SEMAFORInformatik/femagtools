agp_el = 360
local fml = require('fml')
P0 = fml.Point:Create(0,0)
Ls1 = fml.Line:Create(P0,0)
Ls2 = fml.Line:Parallel(Ls1,-${model['slot_width']*1e3/2})
Ls3 = fml.Line:Create(P0,180/m.tot_num_slot)
Ls4 = fml.Line:Create(P0,360/m.tot_num_slot)

Cs1 = fml.Circle:Create(P0,da1/2)
Cs2 = fml.Circle:Create(P0,da1/2+${model['slot_height']*1e3})
Cs3 = fml.Circle:Create(P0,dy1/2)
Cs4 = fml.Circle:Create(P0,da1/2-m.airgap/2)

Ps11 = fml.Point:Intersection(Ls1,Cs1,2)
Ps12 = fml.Point:Intersection(Ls1,Cs2,2)
Ps13 = fml.Point:Intersection(Ls1,Cs3,2)

Ps33 = fml.Point:Intersection(Ls3,Cs3,2)
Ps32 = fml.Point:Intersection(Ls3,Cs2,2)
Ps31 = fml.Point:Intersection(Ls3,Cs1,2)

Ps21 = fml.Point:Intersection(Ls2,Cs1,2)
Ps22 = fml.Point:Intersection(Ls2,Cs2,2)
Ps23 = fml.Point:Intersection(Ls2,Cs3,2)

Ps14 = fml.Point:Intersection(Ls1,Cs4,2)
Ps34 = fml.Point:Intersection(Ls3,Cs4,2)

Ps43 = fml.Point:Intersection(Ls4,Cs3,2)
Ps44 = fml.Point:Intersection(Ls4,Cs4,2)

nc_circle(Ps14.x,Ps14.y, Ps34.x,Ps34.y, agp_el/(2*m.tot_num_slot)+1)

nc_circle(Ps11.x,Ps11.y, Ps21.x,Ps21.y, 0)
nc_circle(Ps21.x,Ps21.y, Ps31.x,Ps31.y, 0)

nc_circle(Ps22.x,Ps22.y, Ps32.x,Ps32.y, 0)

nc_circle(Ps13.x,Ps13.y, Ps23.x,Ps23.y, 0)
nc_circle(Ps23.x,Ps23.y, Ps33.x,Ps33.y, 0)

nc_line(Ps11.x,Ps11.y, Ps12.x,Ps12.y, 0)
nc_line_cont(Ps13.x,Ps13.y, 0)

nc_line(Ps21.x,Ps21.y, Ps22.x,Ps22.y, 0)

nc_line(Ps32.x,Ps32.y, Ps33.x,Ps33.y, 0)

nc_line(Ps11.x,Ps11.y, Ps14.x,Ps14.y, 0)
nc_line(Ps31.x,Ps31.y, Ps34.x,Ps34.y, 0)

create_mesh()

x,y = (Ps11.x+Ps22.x)/2, (Ps11.y+Ps22.y)/2
def_new_sreg(x,y, 'stfe', 'blue')

mirror_nodechains(Ps33.x,Ps33.y, Ps34.x,Ps34.y)
create_mesh()

rotate_copy_nodechains(Ps14.x,Ps14.y, Ps13.x,Ps13.y,
         Ps43.x,Ps43.y, Ps44.x,Ps44.y, m.num_slots-1)

if mcvkey_yoke ~= nil then
  def_mat_fm_nlin(x,y, darkblue, mcvkey_yoke, ${model['rlength']*100})
else
  def_mat_fm(x,y, blue, 500, ${model['rlength']*100})
end

--
x,y = (Ps31.x+Ps32.x)/2,(Ps31.y+Ps32.y)/2
m.wdg_location = -1
m.xcoil_1         = x              --   center coordinate of 1. coil side [mm]
m.ycoil_1         = y              --   center coordinate of 1. coil side [mm]
m.xcoil_2         =          0.000 --   center coordinate of 2. coil side [mm]
m.ycoil_2         =          0.000 --   center coordinate of 2. coil side [mm]
