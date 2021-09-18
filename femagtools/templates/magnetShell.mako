

m.magn_r1    = da1/2  --   Inner Radius Magnet    R1        [mm]   
m.magn_r2    = dy1/2 - ${model.get('yoke_height')*1e3} --   Outer Radius Magnet    R2        [mm]   
m.magn_r4    =       56.5     --   Eccentr. Radius Magnet R4        [mm]   
m.eccentr    =       0.00     --   Eccentricity           E1        [mm]   
m.magn_r3    =       0.00     --   Radius Magnet-end      R3        [mm]   
m.magn_len   = ${model.get('magn_len')*1e2}     --   Magnet Length  [%]   
m.magn_b1    =       100.0     --   Width Magnet outside   B1        [mm]   
m.magn_h1    =       40.0     --   Height Magnet          H1        [mm]   
m.magn_h2    =       14.5     --   Height Magnet          H2        [mm]   
m.magn_rp    =       56.5     --   Radius Housing         Rp        [mm]   
m.magn_a4    =       120.0     --   Angle of R4            A4      [Grad]   
m.magn_rem   =       0.39     --   Remanence              Br         [T]   
m.yoke_height =       6.00     --   Stator Housing height  D         [mm]   
m.magn_ori    =       1.00     --   Orientation: par:1; Pol:2; cos:3        
--m.airgap    =       -1.0     --   Mesh height (2/3 airgap)         [mm]   

m.mcvkey_yoke = mcvkey_yoke                                                    
m.tot_num_sl      =  m.num_poles
m.num_sl_gen      =  m.npols_gen

m.zeroangl  = ${model.get('zeroangl', 0)} 

pre_models("MAGNET_SHELL")
