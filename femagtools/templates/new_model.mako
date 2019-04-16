
  exit_on_error=${model.get('exit_on_error', 'true')}
  exit_on_end=${model.get('exit_on_end')}
  verbosity=${model.get('verbosity', 2)}

model = '${model.get('name')}'
new_model_force(model, '')


<%include file="basic_modpar.mako" />


m.airgap          =     2*ag/3
m.nodedist        =     ${model.stator.get('nodedist',1)}

% if model.get('coord_system', 0) > 0:
m.cood_system     = ${model.get(['coord_system'])} 
% endif
