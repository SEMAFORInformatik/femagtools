-- asyn motor simulation
--
T = ${model.get('bar_temp',20)-20}  --Temperature rise of rotor bar
sigma1, sigma2 = get_dev_data( "cond_conduct" )
tcoeff = 3.9e-3
get_sreg_keys("srkeys") -- get all valid subregions
for j=1, #srkeys do
  srname = get_sreg_data('name', srkeys[j])
  if srname == 'Bar ' then
    sekeys=get_sreg_data("sekeys",srkeys[j])
    ndkeys=get_spel_data("ndkeys",sekeys[1])
    xm, ym = 0, 0
    for i=1, #ndkeys do
      x, y = get_node_data("xy", ndkeys[i])
      xm = xm + x
      ym = ym + y
    end
    x,y=xm/#ndkeys, ym/#ndkeys
    set_mat_cond(x,y, sigma2/(1+tcoeff*T), 100.0 )
  end
end

 post_models("end_wind_leak","leak")
 L_end = math.sqrt(leak[2]^2+leak[3]^2)

-- stator winding resistance per phase
a = ${model.get('num_par_wdgs',1)}  -- parallel winding groups
rl1,rl2 = get_dev_data( "rel_cond_length" )
--[[ TODO: fix m.dia_wire != 0
A_wire = m.dia_wire^2/4*math.pi -- A_ns*m.cufilfactor/m.num_wires
 T = ${model.get('wind_temp',20)-20}  --Temperature rise of stator winding
tcoeff = 3.93e-3 -- temp coeff 1/K
sigma = sigma1/(1+tcoeff*T) -- conductivity 
R_s = Q1/m.num_phases*m.num_wires*(rl1/100)*(m.arm_length/1000)/(sigma*A_wire/1.0e6)/a
--]]
R_s = 0

-- effective rotor bar length (including ring segment)
% if model.get('bar_len', 0):
p = m.num_poles/2
m.num_phases = 3 -- TODO: fix this
length_eff = ${model.get('bar_len')*1e3}
% else:
Dr = da2-m.slot_height
length_eff = rl2/100*m.arm_length+math.pi*Dr/Q2/math.sin(math.pi*p/Q2)
%endif
m.stator_volt     = ${model.get('u1')}       --   Stator windgs (Ph) voltage (RMS) [V]    
m.connect         = ${model.get('wdgcon',0)} --   Wdgs-connect: 0=open;1=star;2=delta     
m.frequency       = ${model.get('f1')}      --   Nominal Stator frequency [Hz]           
m.re_winding      = R_s            --   Stator phase winding resistamce [Ohm]   
m.l_endwindg      = 2*math.pi*m.frequency*L_end          --   Stator ph. end-winding reactance[Ohm]   
m.eff_arm_len     = length_eff     --   Effect. rotor bar length (+endr) [mm]   
m.nphases         = m.num_phases          --   Number of Phases (poc.file)   (>= 2)    
m.num_par_wdgs    = a  --   Number of parallel windings   (>= 0)    
m.slip1           = 0.0            --   Slip 1 [%] (If s1=s2=s3 : s1 )          
m.slip2           = (1-p*${model.get('speed',0)}/m.frequency)*100 --   Slip 2 [%]                              
m.slip3           = 2*m.slip2           --   Slip 3 [%]                              
m.move_action     =          0.000 --   Rot-Motor:0.0; Lin-Motor:2xTaupol[mm]   
m.fc_radius1      = m.fc_radius     --   position [mm] of move path in air gap   

m.pocfilename    = model..'_'..m.num_poles..'p.poc'                                            
run_models("asyn_motor")
