-- Short Circuit
m.re_winding     =     ${model.get('r1',0)}     --   Resistance stator winding  Ra   [Ohm]   
m.l_daxis        =     ${model.get('ld',0)}     --   Inductance                 Ld  [H/mm]   
m.l_qaxis        =     ${model.get('lq',0)}     --   Inductance                 Lq  [H/mm]   
m.l_endwindg     =     ${model.get('l_endwinding',0)}     --   Stator ewdg inductance     Le     [H]   
m.l_external     =     ${model.get('l_external',0)}     --   Stator external inductance Lex    [H]   
m.magn_flux      =     ${model.get('psim',0)} --   Magn. flux (RMS) = Up/omega   [Vs/mm]   
m.arm_length     =     ${model.get('lfe',0)*1e3}     --   Effect. armature length    lm    [mm]   
m.current        =     ${model.get('current',0)}     --   Current (operat. limit) (RMS)     [A]   
m.angl_i_up      =     ${model.get('angl_i_up',0)}     --   Angle current vs. voltage Vp    [Deg]   
m.speed          =     ${model.get('speed',0)*60} --   Speed                         [1/min]   
m.num_pol_pair   =     ${model.get('num_pol_pair',0)}     --   Number of Pole pairs     p  (>= 1)      
m.typ            =     ${model.get('sc_type',3)}     --   Type of sc : 3ph = 3                    
m.initial        =     ${model.get('initial',2)}     --   Initial conditions: noload:1; load:2    
m.simultime      =     ${model.get('simultime',0.1)}     --   Simulation time                   [s]   
m.fc_radius      =     ${model.get('fc_radius',0)*1e3}     --   Radius air-gap center (torq.)    [mm]   
m.num_par_wdgs   =     ${model.get('num_par_wdgs',0)}     --   Number of parallel Windings(>= 1) a     
m.allow_demagn   =     ${model.get('allow_demagn',0)}     --   Allow Demagnetisation:= 1:yes;= 0:no    
m.sim_demagn     =     ${model.get('sim_demagn',0)}     --   Simulate Demagnetisation:1:yes; 0:no    
 
 m.pocfilename   = '${model.get('pocfilename','sin.poc')}'                                  
 
 run_models("shortcircuit")                                            
