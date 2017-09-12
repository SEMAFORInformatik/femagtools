%if model.get('mcvkey',0):
 m.mcvkey_magnet = '${model.get('mcvkey')}'
 m.orient     = ${model.get('orient', 'mpolaniso')}
 m.rlen       = ${model.get('rlen', 1.0)*100}   
 m.magncond   =  ${model.get('magncond', 625000.0)}

%else:
--  Magnet-data
 m.remanenc       =    ${model.get('remanenc')}
 m.relperm        =    ${model.get('relperm')}
 m.spmaweight     =    ${model.get('spmaweight', 7.5e3)*1e-3}
 m.temcoefbr      =    ${model.get('temcoefbr', -0.001)*100}
 m.temcoefhc      =    ${model.get('temcoefhc', -0.001)*100}
 m.magntemp       =    ${model.get('magntemp', 20)}
 m.magncond       =    ${model.get('magncond', 625000.0)}
 m.magnwidth      =    ${model.get('magnwidth', 0.0)*1e3}
 m.magnlength     =    ${model.get('magnlength', 0.0)*1e3}
 
 pre_models("Magnet-data")
 %endif
