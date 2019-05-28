% if model.get('colorgrad_babs',None):
  bmin = ${model.get('bmin', 0.0)}
  bmax = ${model.get('bmax', 0.0)}
  color_gradation( 0,0, 'tot', 'Babs', bmin, bmax, '${model['colorgrad_babs']}')
% endif
% if model.get('colorgrad_demag',None):
  hmin = ${model.get('hmin', 0.0)}
  hmax = ${model.get('hmax', 0.0)}
  color_gradation( 0,0, 'tot', 'demag', hmin, hmax, '${model['colorgrad_demag']}')
% endif

