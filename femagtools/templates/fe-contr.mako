m.hc_min          =         95.000 --   Limit demagnetisa > 0:[%]Hc,<0:[kA/m]   
m.con_hdcopy      =          0.000 --   Hc-copy:Name:auto:0,intact:1, none:-1   
m.b_max           =          2.400 --   Max Induction [T] in colorgradation     
m.b_min           =    ${model.get('move_inside')} --   Move inside: 0 , Move outside: > 0      
m.calc_fe_loss    =          1.000 --   Calc. FE-Loss:0:no, 1:yes, 2:m-output   
m.eval_force      =          0.000 --   Eval. force density > 0, no <= 0        
m.allow_draw      =          1.000 --   Draw Graphics :> 0: yes, 0:  no         
m.fline_dens      =          4.000 --   F-Lines: 1: small, 2: medium, 3:thick   
m.num_flines      =         20.000 --   Number of Field-lines:      < 100 > 2   
m.name_bch_log    =          0.000 --   Name bch-file in Logfile:> 0:yes,0:no   
m.st_size_move    =          0.000 --   Step size move: r/ph:[degr], x/y:[mm]   
m.num_nonl_it     =        300.000 --   Number of nonlinear Iterations   < 99   
m.perm_mode       =          0.000 --   Permeability mode:>0:restore,0:actual   
m.error_perm      =          0.005 --   Rel. Permeability error < 0.1     [%]   
m.allow_demagn    =          0.000 --   Allow Demagnetisation:= 1:yes,= 0:no    
m.maenergy        =          0.000 --   Force from magn energy 1 :yes,= 0:no    
m.el_order_ag     =          1.000 --   El. order in air gap: lin=1: quadr=2    
m.export_scrpt    =          0.000 --   Export parameters in script: yes > 0    

pre_models("FE-contr-data")
