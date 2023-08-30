--  Magnet-data
%if model.get('remanenc',0):
 m.remanenc       =    ${model.get('remanenc')}
%endif
%if model.get('relperm',0):
 m.relperm        =    ${model.get('relperm')}
%endif
%if model.get('spmaweigth',0):
 m.spmaweight     =    ${model.get('spmaweight')*1e-3}
%endif
%if model.get('temcoefbr',0):
 m.temcoefbr      =    ${model.get('temcoefbr')*100}
%endif
%if model.get('temcoefhc',0):
 m.temcoefhc      =    ${model.get('temcoefhc')*100}
%endif
%if model.get('magntemp',0):
 m.magntemp       =    ${model.get('magntemp')}
%endif
%if model.get('magncond',0):
 m.magncond       =    ${model.get('magncond')}
%endif
%if model.get('rlen',0):
 m.rlen           =    ${model.get('rlen')*100}
%endif
%if model.get('mcvkey',0):
 m.mcvkey_magnet = '${model.get('mcvkey')}'
%if model.get('orient', 0):
 m.orient     = ${model.get('orient', 'm.cartiso')}
%endif
%else:
%if model.get('magnwidth',0):
 m.magsegwid      =    ${model.get('magnwidth')*1e3}
%elif model.get('magnsegwidth',0):
 m.magsegwid      =    ${model.get('magnsegwidth')*1e3}
%endif
%if model.get('magnlength',0):
 m.magseglen     =    ${model.get('magnlength')*1e3}
%endif
%if model.get('magnseglength',0):
 m.magseglen     =    ${model.get('magnseglength')*1e3}
%endif
%endif

pre_models("Magnet-data")
