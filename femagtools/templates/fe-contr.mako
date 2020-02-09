m.hc_min          =  ${'%12.3f' % model.get('hc_min', 95.0)} --   Limit demagnetisa > 0:[%]Hc,<0:[kA/m]   
m.con_hdcopy      =  ${'%12.3f' % model.get('con_hdcopy', 0)} --   Hc-copy:Name:auto:0,intact:1, none:-1   
m.b_max           =  ${'%12.3f' % model.get('b_max', 2.4)} --   Max Induction [T] in colorgradation     
m.b_min           =  ${'%12.3f' % model.get('move_inside')} --   Move inside: 0 , Move outside: > 0      
m.calc_fe_loss    =  ${'%12.3f' % model.get('calc_fe_loss', 1)} --   Calc. FE-Loss:0:no, 1:yes, 2:m-output   
m.eval_force      =  ${'%12.3f' % model.get('eval_force', 0)} --   Eval. force density > 0, no <= 0        
m.allow_draw      =  ${'%12.3f' % model.get('allow_draw', 1)} --   Draw Graphics :> 0: yes, 0:  no         
m.fline_dens      =  ${'%12.3f' % model.get('fline_dens', 4)} --   F-Lines: 1: small, 2: medium, 3:thick   
m.num_flines      =  ${'%12.3f' % model.get('num_flines', 20)} --   Number of Field-lines:      < 100 > 2   
m.name_bch_log    =  ${'%12.3f' % model.get('name_bch_log', 0)} --   Name bch-file in Logfile:> 0:yes,0:no   
m.st_size_move    =  ${'%12.3f' % model.get('st_size_move', 0)} --   Step size move: r/ph:[degr], x/y:[mm]   
m.num_nonl_it     =  ${'%12.3f' % model.get('num_nonl_it', 300)} --   Number of nonlinear Iterations   < 99   
m.perm_mode       =  ${'%12.3f' % model.get('perm_mode', 0)} --   Permeability mode:>0:restore,0:actual   
m.error_perm      =  ${'%12.3f' % model.get('error_perm', 0.005)} --   Rel. Permeability error < 0.1     [%]   
m.allow_demagn    =  ${'%12.3f' % model.get('allow_demagn', 0)} --   Allow Demagnetisation:= 1:yes,= 0:no    
m.maenergy        =  ${'%12.3f' % model.get('maenergy', 0)} --   Force from magn energy 1 :yes,= 0:no    
m.el_order_ag     =  ${'%12.3f' % model.get('el_order_ag', 1)} --   El. order in air gap: lin=1: quadr=2    
m.export_scrpt    =  ${'%12.3f' % model.get('export_scrpt', 0)} --   Export parameters in script: yes > 0    

pre_models("FE-contr-data")
