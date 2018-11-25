
--  FE-Losses
m.basfreq         =      ${model.get('basfreq', 50.0)}
m.basind          =      ${model.get('basind', 1.5)}
m.ch              =      ${model.get('ch', 4.0)}
m.cw              =      ${model.get('cw', 2.0)}
m.hyscoef         =      ${model.get('hyscoef', 1.0)}
m.edycoef         =      ${model.get('edycoef', 2.0)}
m.indcoef         =      ${model.get('indcoef', 2.0)}
m.ffactor         =      ${model.get('ffactor', 1.0)}
m.spweight        =      ${model.get('spweight', 7.65)}
m.fillfact        =      ${model.get('fillfact', 1.0)}
m.emodul          =      ${model.get('emodul', 210e9)}
m.poison          =      ${model.get('poison', 0.3)}
m.dampfact        =      ${model.get('dampfact', 0.0)}
m.thcond          =      ${model.get('thcond', 30.0)}
m.thcap           =      ${model.get('thcap', 480.0)}

 pre_models("FE-Losses-1")
 pre_models("FE-Losses-2")
