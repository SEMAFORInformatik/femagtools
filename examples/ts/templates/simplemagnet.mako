Lr1 = fml.Line:Create(P0,0)
Lr3 = fml.Line:Create(P0,180/m.num_poles)
Lr2 = fml.Line:Parallel(Lr3,${model['magnet_width']*1e3/2})
Lr4 = fml.Line:Create(P0,360/m.num_poles)

Cr1 = fml.Circle:Create(P0,da2/2)
Cr2 = fml.Circle:Create(P0,da2/2-${model['magnet_height']*1e3})
Cr3 = fml.Circle:Create(P0,dy2/2)
Cr4 = fml.Circle:Create(P0,da2/2+m.airgap/2)

Pr11 = fml.Point:Intersection(Lr1,Cr1,2)
Pr12 = fml.Point:Intersection(Lr1,Cr2,2)
Pr13 = fml.Point:Intersection(Lr1,Cr3,2)

Pr21 = fml.Point:Intersection(Lr2,Cr1,2)
Pr22 = fml.Point:Intersection(Lr2,Cr2,2)

Pr31 = fml.Point:Intersection(Lr3,Cr1,2)
Pr32 = fml.Point:Intersection(Lr3,Cr2,2)
Pr33 = fml.Point:Intersection(Lr3,Cr3,2)

Pr14 = fml.Point:Intersection(Lr1,Cr4,2)
Pr34 = fml.Point:Intersection(Lr3,Cr4,2)
Pr44 = fml.Point:Intersection(Lr4,Cr4,2)

nc_circle(Pr11.x,Pr11.y, Pr21.x,Pr21.y, 0)
nc_circle(Pr21.x,Pr21.y, Pr31.x,Pr31.y, 0)

nc_circle(Pr12.x,Pr12.y, Pr22.x,Pr22.y, 0)
nc_circle(Pr22.x,Pr22.y, Pr32.x,Pr32.y, 0)

nc_circle(Pr13.x,Pr13.y, Pr33.x,Pr33.y, 0)

nc_circle(Pr14.x,Pr14.y, Pr34.x,Pr34.y, agp_el/(2*m.num_poles)+1)

nc_line(P0.x,P0.y, Pr13.x,Pr13.y, 0)
nc_line(Pr13.x,Pr13.y, Pr12.x,Pr12.y, 0)
nc_line(Pr12.x,Pr12.y, Pr11.x,Pr11.y, 0)

nc_line(P0.x,P0.y, Pr33.x,Pr33.y, 0)
nc_line(Pr33.x,Pr33.y, Pr32.x,Pr32.y, 0)

nc_line(Pr22.x,Pr22.y, Pr21.x,Pr21.y, 0)

nc_line(Pr11.x,Pr11.y, Pr14.x,Pr14.y, 0)
nc_line(Pr31.x,Pr31.y, Pr34.x,Pr34.y, 0)

create_mesh()

x,y = (Pr12.x+Pr32.x)/2, (Pr12.y+Pr32.y)/2
fml.def_new_sreg(x,y, 'rofe', blue)

x,y = (Pr13.x+Pr33.x)/2, (Pr13.y+Pr33.y)/2
fml.def_new_sreg(x,y, 'wefe', cyan)

mirror_nodechains(Pr34.x,Pr34.y, P0.x,P0.y)

create_mesh()

rotate_copy_nodechains(P0.x,P0.y, Pr14.x,Pr14.y,
              Pr44.x,Pr44.y, P0.x,P0.y, m.npols_gen-1)

x,y = (Pr12.x+Pr32.x)/2,(Pr12.y+Pr32.y)/2
if mcvkey_yoke ~= nil then
  def_mat_fm_nlin(x,y, darkblue, mcvkey_yoke, 100)
else
  def_mat_fm(x,y, blue, 500, 100)
end

x,y = (Pr13.x+Pr33.x)/2,(Pr13.y+Pr33.y)/2
def_mat_fm(x,y, blue, 200, 100)

for i = 1,m.npols_gen,2 do
  Br = 1.2
  muer = 1.0
  a = (i-0.5)*360/m.num_poles
  x,y = pd2c((da2-${model['magnet_height']*1e3})/2,a)
  def_mat_pm(x,y, 'red', Br, muer, a, m.parallel, 0.0, ${model['rlength']*100})
  a = (i+0.5)*360/m.num_poles
  x,y = pd2c((da2-${model['magnet_height']*1e3})/2,a)
  def_mat_pm(x,y, 'green', Br, muer, 180+a, m.parallel, 0.0, ${model['rlength']*100})
end
