-- calculate static or dynamic eccentricity
-- NOTE: requires full model
--
static = 0
dynamic = 1
m.disp_info   =  ${model.get('type', 'static')}  -- Displace: Stator: 0, Rotor: 1

da1 = ${model.get('bore_diam',0)*1e3}
ag = ${model.get('airgap',0)*1e3}

if m.b_min == 0 then -- move inside: outer stator
  if m.disp_info == static then
    m.fc_radius = da1/2 - ag/6
  else
    m.fc_radius = da1/2 - 5*ag/6
  end
else
  if m.disp_info == static then
    m.fc_radius = da1/2 + ag/6
  else
    m.fc_radius = da1/2 + 5*ag/6
  end
end
m.disp_radius =  ${model.get('ecc', 0)*1e3}  -- Displacement: r [mm]
m.disp_phi    =  0.0         -- Displacement: phi [Degr]

pre_models("Displ_Stat/Rot")
-- must correct fc_radius
if m.b_min == 0 then -- move inside: outer stator
  m.fc_radius = da1/2 - ag/2
else
  m.fc_radius = da1/2 + ag/2
end