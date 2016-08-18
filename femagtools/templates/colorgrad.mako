% if model.get('colorgrad_babs',None):
  color_gradation( 0,0, tot, Babs, 0, 2.4, '${model['colorgrad_babs']}')
% endif
% if model.get('colorgrad_demag',None):
  color_gradation( 0,0, tot, demag, 0, 2.4, '${model['colorgrad_demag']}')
% endif

