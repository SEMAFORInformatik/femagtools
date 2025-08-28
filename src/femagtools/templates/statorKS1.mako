
m.yoke_diam   = dy1
m.inside_diam = da1

m.wdg_location = -1 -- for gen_windings

m.slot_angle = ${model['slot_angle']*1e3}
m.slot_height = ${model['slot_height']*1e3}
m.slot_h1     = ${model['slot_h1']*1e3}
m.slot_h2     = ${model['slot_h2']*1e3}
m.slot_bk     = ${model['slot_bk']*1e3}
m.slot_width  = ${model['slot_width']*1e3}
m.slot_r1     = ${model['slot_r1']*1e3}
m.slot_r2     = ${model['slot_r2']*1e3}
m.middle_line = ${model.get('middle_line',0)}
m.tooth_width = ${model.get('tooth_width',0)*1e3}
m.slot_topwidth = ${model.get('slot_topwidth',0)*1e3}

m.zeroangl    = ${model.get('zeroangle',0)}
m.rlength     = ${model.get('rlength',1)*100}

% if model.get('ac_loss', False):
  m.ac_loss = ${model.get('ac_loss', 6)}
% endif
m.mcvkey_yoke = mcvkey_yoke

pre_models("STATOR_KS1")

if mcvkey_teeth ~= nil then
  if m.inside_diam > m.yoke_diam then
     r = (m.inside_diam - m.slot_height)/2
  else
     r = (m.inside_diam + m.slot_height)/2
  end
  x0, y0 = pr2c(r, 2*math.pi/m.tot_num_slot + m.zeroangl/180*math.pi)
   def_mat_fm_nlin(x0, y0, "blue", mcvkey_teeth, m.rlength)
end

%if model.get('thcond', 0) and model.get('thcap', 0):
stator_thcond = ${model.get('thcond', 24)}
stator_thcap = ${model.get('thcap', 480)}
stator_density = ${model.get('density', 7700)}
%endif
