

m.magn_r1    = da1/2  --   Inner Radius Magnet    R1 [mm]
m.magn_r2    = dy1/2 - ${model.get('yoke_height')*1e3} -- Outer Radius Magnet R2 [mm]
m.magn_r4    = ${model.get('magn_r4')*1e3}  -- Eccentr. Radius Magnet R4 [mm]
m.eccentr    = ${model.get('eccentr')*1e3}  -- Eccentricity  E1 [mm]   
m.magn_r3    = ${model.get('eccentr')*1e3}  -- Radius Magnet-end R3  [mm]
m.magn_len   = ${model.get('magn_len')*1e2} -- Magnet Length  [%]   
m.magn_b1    = ${model.get('magn_b1')*1e3}  -- Width Magnet outside B1 [mm]
m.magn_h1    = ${model.get('magn_h1')*1e3}  -- Height Magnet H1 [mm]
m.magn_h2    = ${model.get('magn_h2')*1e3}  -- Height Magnet H2 [mm]
m.magn_rp    = ${model.get('magn_rp')*1e3}  -- Radius Housing         Rp [mm]
m.magn_a4    = ${model.get('magn_a4')}  --   Angle of R4 A4 [Grad]   
m.magn_rem   =       0.0     --   (ignored) Remanence  Br [T]   
m.yoke_height = ${model.get('yoke_height')*1e3} --   Stator Housing height D [mm]
m.magn_ori    = ${model.get('magn_ori')} --  Orientation: par:1; Pol:2; cos:3
--m.airgap    =  -1.0  --   Mesh height (2/3 airgap) [mm]

m.mcvkey_yoke = mcvkey_yoke                                                    
m.tot_num_sl  =  m.num_poles
m.num_sl_gen  =  m.npols_gen

m.zeroangl  = ${model.get('zeroangl', 0)} 

pre_models("MAGNET_SHELL")
