
m.stator_diam     =      da1
m.inside_diam     =      dy2
m.airgap          =     ag
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
Q1 = m.tot_num_sl
Q2 = ${model['num_slots']}
m.tot_num_sl      = Q2
m.num_sl_gen      = Q2 * m.npols_gen/m.num_poles
                                  --    Number of teeth be generated
                                  
m.slot_height = m.slot_h32 + m.slot_h42 + m.slot_h52 +
               m.slot_h62 + m.slot_h72
m.zeroangl        = ${model.get('zeroangl',0)}
m.nodedist        =   ${model.get('nodedist',1)}

pre_models("ROTOR_ASYN")                                       

-- yoke material

if mcvkey_yoke ~= 'dummy' then
   m.rlength         =     ${model.get('rlength', 1)*100}  
   x0, y0 = pd2c(dy2/2+0.1, 360/m.num_sl_gen/2+m.zeroangl)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_yoke, m.rlength)
end

-- rename subregion to Bar
r=da2/2 - m.slot_hs2 - m.slot_height/2
phi=math.pi/Q2+m.zeroangl
x,y=pr2c(r,phi)
delete_sreg(x, y)
x,y=pr2c(r,0.99*phi)
def_new_subreg( x,y, 'Bar', violet )
x,y = pr2c(r,1.01*phi)	
add_to_subreg( x, y )
for i=2, m.num_sl_gen do
  phi=(i-0.5)*360/Q2 + m.zeroangl
  for j = -1,1,2 do
    x,y = pd2c(r,phi+j*0.01)	
    add_to_subreg( x, y )
  end
 end
