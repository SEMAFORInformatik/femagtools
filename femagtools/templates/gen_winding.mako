% if model.windings['wdgtype'] == 'CMM':
--  GEN_CMM_WDG

m.s_or_w_windg    = ${model.windings.get('s_or_w_windg', 1)} --   Lap/wave/frog-leg winding = 1/2/3
m.numpl_sw        = ${model.windings.get('numpl_sw', 1)} --   Number of plex of lap winding: m_l
m.numpl_ww        = ${model.windings.get('numpl_ww', 1)} --   Number of plex of wave winding: m_w
m.foreward        = ${model.windings.get('foreward', 1)} --   Progressive = 1 / retrogressive = 2
m.num_wires       = ${model.windings.get('num_wires', 1)} --   Number of wires per slot side w_sp
m.current         =       0     --   Armat.-Wdg-Current[A] or flux[Vs/mm]
m.wind_type       =       1     --   Wdg-coil:1=w&cur;2=w&flux;3=bar&cur
m.num_layers      = ${model.windings.get('num_layers', 1)} --   Number of coil sides per slot layer:u
m.pitch_fact      = ${model.windings.get('pitch_fact', 1)} --   Short pitch factor : beta_v
m.dc_ac           = ${model.windings.get('dc_ac', 0)}     --   Current: DC: 0; AC: 1

 pre_models("GEN_CMM_WDG")
-- gen cmm
pre_models("GEN_CMM"); -- autm. generate cmm File (requires femag rel >=9.2)"
% elif 'wdgfile' in model.windings:
def_new_wdg('${model.windings.get("wdgfile")}')
pre_models("gen_pocfile")
% else:
--  Gen_winding
if m.xcoil_1 ~= nil then
  m.wdg_location = 1 --stator
end

m.num_phases      =  ${model.get(['windings','num_phases'])}
m.num_layers      =  ${model.get(['windings','num_layers'])}
m.num_wires       =  ${model.get(['windings','num_wires'])}
m.coil_span       =  ${model.get(['windings','coil_span'])}
% if 'num_poles' in model.windings:
m.num_poles       =  ${model.get(['windings','num_poles'])}
% endif
% if model.get('move_action', 0) == 0:
m.current         =   0.0
m.mat_type        =   1.0 -- rotating
m.wind_type       =   1.0 -- winding & current
m.win_asym        =   1.0 -- sym

m.curr_inp        =   0.0 -- const
m.dq_offset       =   0

% if model.get(['stator', 'num_slots_gen']) == 1:
def_new_wdg(m.xcoil_1, m.ycoil_1, "green", "1", m.num_wires, 0.0, "wi")
add_to_wdg(m.xcoil_2, m.ycoil_2, "wsamekey", "wo", "wser")
% else:
pre_models("Gen_winding")
pre_models("gen_pocfile")
% endif
% else:  # move_action > 0
color={"green", "yellow", "magenta", "lightgrey", "darkred", "skyblue", "violet"}
wkey={0,0,0,0,0,0}

bz = m.width_bz
sw = m.slot_width
ys = m.slot_height/2
wdgscheme = ${model.windings.get('wdgscheme', '{}')}
-- TODO: m.middle_line = 1 only
for l=1, #wdgscheme do
  for z=1, #wdgscheme[l] do
    for i=1, #wdgscheme[l][z] do
      k = wdgscheme[l][z][i]
      if math.abs(k) < m.num_sl_gen+1 then
        xs = (2*math.abs(k)-1)*bz/2
        dir = 1
        if k < 0 then
          dir = -1
        end
        if wkey[z] == 0 then
         wdg = "wdg"..z
         wkey[z]=def_new_wdg(xs + sw/4, ys, color[z], wdg, m.num_wires, 0, dir)
        else
          if l == 1 then
             add_to_wdg(xs + sw/4, ys, wkey[z], dir, "wser")
          else
             add_to_wdg(xs - sw/4, ys, wkey[z], dir, "wser")
          end
        end
      end
    end
  end
end

----------------------------------------------------------------------
-- create poc-File
----------------------------------------------------------------------
base_angle = 0.0
period = 2*m.tot_num_slot*m.width_bz/m.num_poles
f = assert(io.open(model..'_'..m.num_poles.."p.poc","w"));
  f:write(string.format("          %i\n",m.num_phases));
  for i = 1, m.num_phases do
    if i<10      then f:write(string.format("          %i\n",i));
    elseif i<100 then f:write(string.format("         %i\n",i));
    end
  end
  for i = 0,m.num_phases-1 do
    angle = base_angle+i*360.0/m.num_phases
    if     angle<10   then f:write(string.format("   %10.8f\n",angle));
    elseif angle<100  then f:write(string.format("   %10.7f\n",angle));
    elseif angle<1000 then f:write(string.format("   %10.6f\n",angle));
    end
  end
  if     period<10   then f:write(string.format("   %10.8f\n",period));
  elseif period<100  then f:write(string.format("   %10.7f\n",period));
  elseif period<1000 then f:write(string.format("   %10.6f\n",period));
  end
  f:write("sin\n");
  f:write("   0.0\n");
  f:write("          0\n");
io.close(f);

%endif
% if 'num_poles' in model.windings:
m.num_poles       =  ${model.poles}
% endif
% endif

% if 'thcap' in model.windings:
-- Thermal Material properties
if m.slot_height ~= nil then
  -- FEMAG slot model
  -- TODO: slot model from user
  rw = da1/2 +(m.slot_height-m.slot_h1)/2
  dw = 0
  dr = 0
  if m.middle_line == 1 then
    dw = 1/60
  elseif m.middle_line == 2 then
    dr = 1
  end
  lamCu = 400
  capCu = 385
  da = 1.0785
  dCu = 1.0
  lam = lamCu*(dCu/(da-dCu)+(da-dCu)/da)
  cap = capCu*da^2/(dCu^2*math.pi/4)
  for i=1,m.num_sl_gen do
    a = (2*i-1)*math.pi/m.tot_num_sl + m.zeroangl/180*math.pi
    xwl,ywl = pr2c(rw+dr,a+dw)
    def_mat_therm(xwl,ywl,'yellow',8920,lam,cap,1)
    if m.middle_line > 0 then
      xwr,ywr = pr2c(rw-dr,a-dw)
      def_mat_therm(xwr,ywr,'yellow',8920,lam,cap,1)
    end
  end
end
%endif