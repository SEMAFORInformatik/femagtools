
m.stator_diam     =      da1
m.inside_diam     =      dy2

m.slot_bs2        =     ${model['slot_bs2']*1e3}
m.slot_hs2        =     ${model['slot_hs2']*1e3}
m.slot_b32        =     ${model['slot_b32']*1e3}
m.slot_h32        =     ${model['slot_h32']*1e3}
m.slot_b42        =     ${model['slot_b42']*1e3}
m.slot_h42        =     ${model['slot_h42']*1e3}
m.slot_b52        =     ${model['slot_b52']*1e3}
m.slot_b62        =     ${model['slot_b62']*1e3}
m.slot_h52        =     ${model['slot_h52']*1e3}
m.slot_h62        =     ${model['slot_h62']*1e3}
m.slot_h72        =     ${model['slot_h72']*1e3}
m.mcvkey_yoke     = mcvkey_yoke
Q1 = m.tot_num_sl
Q2 = ${model['num_slots']}
m.tot_num_sl      = Q2
% if model.get('num_slots_gen', 0):
m.num_sl_gen  =   ${model['num_slots_gen']}
% else:
m.num_sl_gen  =   Q2 * m.npols_gen/m.num_poles
% endif

hs = { m.slot_h32, m.slot_h42,
       m.slot_h52 + m.slot_h62 + m.slot_h72 }
m.zeroangl        = ${model.get('zeroangle',0)}
m.nodedist        =   ${model.get('nodedist',1)}

pre_models("ROTOR_ASYN")

-- yoke material

if mcvkey_yoke ~= 'dummy' then
   m.rlength         =     ${model.get('rlength', 1)*100}
   x0, y0 = pd2c(dy2/2+0.1, 180/m.num_sl_gen+m.zeroangl)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_yoke, m.rlength)
end

-- rename subregion to Bar
dphi=1e-3
r=da2/2 - m.slot_hs2 - hs[1] - hs[2] - hs[3]/2
phi=180/Q2+m.zeroangl
x,y=pd2c(r,phi)
delete_sreg(x, y)
x,y=pd2c(r,phi-dphi)
def_new_subreg( x,y, 'Bar', violet )
x,y = pd2c(r,phi+dphi)
add_to_subreg( x, y )

if(hs[2]>0.1) then
  r=da2/2 - m.slot_hs2 - hs[1] - hs[2]/2
  x,y=pd2c(r,phi-dphi)
  add_to_subreg( x,y )
  x,y = pd2c(r,phi+dphi)
  add_to_subreg( x, y )
end
if(hs[1]>0) then
  r=da2/2 - m.slot_hs2 - hs[1]/2
  x,y=pd2c(r,phi-dphi)
  add_to_subreg( x,y )
  x,y = pd2c(r,phi+dphi)
  add_to_subreg( x, y )
end

for i=2, m.num_sl_gen do
  phi=(i-0.5)*360/Q2 + m.zeroangl
  for j = -1,1,2 do
    r=da2/2 - m.slot_hs2 - hs[1] - hs[2] - hs[3]/2
    x,y = pd2c(r,phi+j*dphi)
    add_to_subreg( x, y )
    if(hs[1]>0.1) then
      r=da2/2 - m.slot_hs2 - hs[1]/2
      x,y = pd2c(r,phi+j*dphi)
      add_to_subreg( x, y )
    end
    if(hs[2]>0.1) then
      r=da2/2 - m.slot_hs2 - hs[1] - hs[2]/2
      x,y=pd2c(r,phi+j*dphi)
      add_to_subreg( x,y )
    end
  end
 end
