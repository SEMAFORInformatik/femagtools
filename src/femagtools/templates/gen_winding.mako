% if model.winding['wdgtype'] == 'CMM':
--  GEN_CMM_WDG

m.s_or_w_windg    = ${model.winding.get('s_or_w_windg', 1)} --   Lap/wave/frog-leg winding = 1/2/3
m.numpl_sw        = ${model.winding.get('numpl_sw', 1)} --   Number of plex of lap winding: m_l
m.numpl_ww        = ${model.winding.get('numpl_ww', 1)} --   Number of plex of wave winding: m_w
m.foreward        = ${model.winding.get('foreward', 1)} --   Progressive = 1 / retrogressive = 2
m.num_wires       = ${model.winding.get('num_wires', 1)} --   Number of wires per slot side w_sp
m.current         =       0     --   Armat.-Wdg-Current[A] or flux[Vs/mm]
m.wind_type       =       1     --   Wdg-coil:1=w&cur;2=w&flux;3=bar&cur
m.num_layers      = ${model.winding.get('num_layers', 1)} --   Number of coil sides per slot layer:u
m.pitch_fact      = ${model.winding.get('pitch_fact', 1)} --   Short pitch factor : beta_v
m.dc_ac           = ${model.winding.get('dc_ac', 0)}     --   Current: DC: 0; AC: 1

 pre_models("GEN_CMM_WDG")
-- gen cmm
pre_models("GEN_CMM"); -- autm. generate cmm File (requires femag rel >=9.2)"
% elif 'wdgfile' in model.winding:
def_new_wdg('${model.winding.get("wdgfile")}')
pre_models("gen_pocfile")
% else:
--  Gen_winding
if m.xcoil_1 ~= nil then
  m.wdg_location = 1 --stator
end

m.num_phases      =  ${model.get(['winding','num_phases'])}
m.num_layers      =  ${model.get(['winding','num_layers'])}
m.num_wires       =  ${model.get(['winding','num_wires'])}
m.coil_span       =  ${model.get(['winding','coil_span'])}
% if 'num_poles' in model.winding:
m.num_poles       =  ${model.get(['winding','num_poles'])}
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
wdgscheme = ${model.winding.get('wdgscheme', '{}')}
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
%endif
% if 'num_poles' in model.winding:
m.num_poles       =  ${model.poles}
% endif
% endif

% if 'thcap' in model.winding:
-- Thermal Material properties

conductor_density = ${model.winding['spmaweight']*1e3}
conductor_thcond = ${model.winding['thcond']}
conductor_thcap = ${model.winding['thcap']}

--[[
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

  density = ${model.winding['spmaweight']*1e3}
  lam = ${model.winding['thcond']}
  cap = ${model.winding['thcap']}

  for i=1,m.num_sl_gen do
    a = (2*i-1)*math.pi/m.tot_num_sl + m.zeroangl/180*math.pi
    xwl,ywl = pr2c(rw+dr,a+dw)
    def_mat_therm(xwl,ywl,'yellow',density,lam,cap,1)
    if m.middle_line > 0 then
      xwr,ywr = pr2c(rw-dr,a-dw)
      def_mat_therm(xwr,ywr,'yellow',density,lam,cap,1)
    end
  end
end
]]--
%endif
