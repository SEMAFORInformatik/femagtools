--  CU-Losses-                           
m.cufilfact       =      ${model['windings'].get('cufilfact', 0.45)}
m.culength        =      100*${model['windings'].get('culength', 1.4)}
m.cuconduct       =      56e6
m.numlayers       =      ${model['windings']['num_layers']}
m.conheight       =      0.0
m.contemp         =      ${model.get('wind_temp',20.0)}
m.emodul          =      210.0 
m.poison          =      0.3
m.dampfact        =      0.0
m.thcond          =      30.0
m.thcap           =      480.0
m.slot_indul      =      ${model['windings'].get('slot_indul',0.0)*1e3}
m.dia_wire        =      0.0
m.num_wire        =      0.0
 
pre_models("CU-Losses-1") -- outside
pre_models("CU-Losses-2") -- inside
