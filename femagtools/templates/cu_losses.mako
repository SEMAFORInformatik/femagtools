
--  CU-Losses-
m.cufilfact       =  ${model.get('cufilfact', 0.45)}
m.culength        =  ${model.get('culength', 1.4)*100} -- rel wire length %
m.cuconduct       =  ${model.get('cuconduct', 56e6)}  -- el. conductivity S/m
m.numlayers       =  ${model.get('num_layers', 1)}
m.conheight       =  ${model.get('conheight', 0.0)}
m.contemp         =  ${model.get('wind_temp',20.0)} -- conductor temperature °C
m.emodul          =  ${model.get('emodul', 210e9)}
m.poison          =  ${model.get('poisson', 0.3)}
m.dampfact        =  ${model.get('dampfact', 0.0)}
m.thcond          =  ${model.get('thcond', 30.)}
m.thcap           =  ${model.get('thcap', 480.0)}
m.slot_indul      =  ${model.get('slot_indul',0.0)*1e3}
m.dia_wire        =  ${model.get('dia_wire',0.0)*1e3} -- wire diameter mm
m.num_wire        =  ${model.get('num_wire',0.0)}
% if model.get('winding_inside', False):
-- conductor properties
m.conduct         = m.cuconduct
m.conrelperm      =  ${model.get('relperm', 100)} -- rel permeability
m.contecoef       =  ${model.get('tempcoef', 0)} -- conduct. temperature coeff 1/K
m.spconweight     =  ${model.get('spmaweight', 7.6)} -- mass density g/cm³
m.relconlength    =  ${model.get('rlen', 1)*100} -- rel cond length %
pre_models('conduct-data')
pre_models("CU-Losses-2") -- inside
% else:
pre_models("CU-Losses-1") -- outside
% endif
