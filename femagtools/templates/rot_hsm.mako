m.inside_diam   = dy2
m.rotor_diam    = da2+2*ag

m.gap_pol_shaft = ${model['gap_pol_shaft']*1e3}
m.core_height   = ${model['core_height']*1e3}
m.pole_height   = ${model['pole_height']*1e3}
m.pole_rad      = ${model['pole_rad']*1e3}
m.core_width2   = ${model['core_width2']*1e3}
m.core_width1   = ${model['core_width1']*1e3}
m.pole_width_r  = ${model['pole_width_r']*1e3}
m.pole_width    = ${model['pole_width']*1e3}
m.slot_width    = ${model['slot_width']*1e3}
m.slot_height   = ${model['slot_height']*1e3}
m.damper_diam   = ${model['damper_diam']*1e3}
m.damper_div    = ${model['damper_div']*1e3}
m.mcvkey_yoke   = mcvkey_yoke

m.zeroangl      = ${model.get('zeroangle', 0)}
m.tot_num_sl    = m.num_poles
m.num_sl_gen    = m.npols_gen
m.airgap  = ag

pre_models("ROT_HSM")

-- Yoke material
if mcvkey_yoke ~= 'dummy' then
  x, y = pd2c( m.inside_diam/2+1, 360/m.num_poles/2)
  rellen = 100 -- relative length in %
  def_mat_fm_nlin(x, y, "blue", mcvkey_yoke, rellen)
end
