--  Magnet-data
 m.remanenc       =    ${model.get('remanenc')}
 m.relperm        =    ${model.get('relperm')}
 m.spmaweight     =    ${model.get('spmaweight')}
 m.temcoefbr      =    ${model.get('temcoefbr',-0.1)*100}
 m.temcoefhc      =    ${model.get('temcoefhc',-0.1)*100}
 m.magntemp       =    ${model.get('magntemp',20)}
 m.magncond       =    ${model.get('magncond')}
 m.magnwidth      =    ${model.get('magnwidth')*1e3}
 m.magnlength     =    ${model.get('magnlength')*1e3}
 
 pre_models("Magnet-data")
