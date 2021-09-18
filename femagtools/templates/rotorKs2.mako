

m.rotor_diam    = da2
m.inside_diam   = dy2
m.slot_angle    = ${model.get('slot_angle')*1e3}  -- Angle ALPHA [degr]  
m.slot_height   = ${model.get('slot_height')*1e3} -- Total height H [mm]
m.slot_topwidth = ${model.get('slot_topwidth')*1e3}  -- Slot width top B [mm]
m.slot_width    = ${model.get('slot_width')*1e3}   -- Slot width SW [mm]
m.slot_h1       = ${model.get('slot_h1')*1e3} -- Distance H1 [mm]
m.slot_h2       = ${model.get('slot_h2')*1e3} -- Distance H2 [mm]   
m.slot_r1       = ${model.get('slot_r1')*1e3} -- Radius R1 [mm]
m.slot_r2       = ${model.get('slot_r2')*1e3} -- Radius R2 [mm]   
m.middle_line   = ${model.get('middle_line')} -- Centerline: no: 0; vertic: 1; horiz:2   
m.nodedist        = 1.00     --   Rel. node distance in iron > 0.5        
m.airgap          = -1     --   Airgap 2/3 height:>0:3; <0:4 rows[mm]   
m.tot_num_sl      =  25   --    Total number of teeth on 2pi            
m.num_sl_gen      =  25   --    Number of teeth be generated            

m.zeroangl        = ${model.get('zeroangl',0)}
m.nodedist        =   ${model.get('nodedist',1)}

m.mcvkey_yoke = mcvkey_yoke                                                    

pre_models("ROTOR_KS2");  

if mcvkey_yoke ~= 'dummy' then
  m.rlength         =     ${model.get('rlength', 1)*100}  
--   x0, y0 = pd2c(dy2/2+0.1, 360/m.num_sl_gen/2+m.zeroangl)
--   def_mat_fm_nlin(x0, y0, "blue", mcvkey_yoke, m.rlength)

  x, y = pd2c(m.rotor_diam/2-m.slot_height/2,0)
  def_mat_fm_nlin(x, y, "blue", mcvkey_yoke, m.rlength)
  
  x, y = pd2c(m.inside_diam/2+(m.rotor_diam/2-m.inside_diam/2-m.slot_height)/2,0)
  def_mat_fm_nlin(x, y, "blue", mcvkey_yoke, m.rlength)

end
