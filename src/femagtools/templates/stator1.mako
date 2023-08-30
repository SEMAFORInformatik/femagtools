m.yoke_rad   = dy1/2
m.inside_rad = da1/2

m.tip_rh2         =     ${model['tip_rh2']*1e3}
m.tip_rh1         =     ${model['tip_rh1']*1e3}
m.slot_rf1        =     ${model['slot_rf1']*1e3}
m.slot_rf2        =     m.slot_rf1
m.yoke_rad2       =     m.yoke_rad
m.slot_width      =     ${model['slot_width']*1e3}
m.tooth_width     =     ${model['tooth_width']*1e3}
m.zeroangl        =     ${model['zeroangle']}
m.rlength         =     ${model['rlength']*100}  

 m.mcvkey_yoke    =   mcvkey_yoke
 m.wdg_location   =   -1.0 -- stator (internal values)

pre_models( "STATOR_1")

if mcvkey_teeth ~= nil then
  x0, y0 = pr2c( (m.tip_rh1 + m.slot_rf1)/2, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_teeth, m.rlength)
end
