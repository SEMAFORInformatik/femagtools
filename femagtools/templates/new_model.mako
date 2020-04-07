
  exit_on_error=${model.get('exit_on_error', 'false')}
  exit_on_end=${model.get('exit_on_end')}
  verbosity=${model.get('verbosity', 2)}

model = '${model.get('name')}'
description = '${model.get('description','')}'
% if model.get('scratch_mode', 0):
new_model(model, '', 'scratch')
% else:
new_model_force(model, '')
%endif

<%include file="basic_modpar.mako" />

% if isinstance(model.get(['bore_diam']), list):
m.airgap          =   2*ag[2]/3
% else:
m.airgap          =     2*ag/3
% endif
m.nodedist        =     ${model.stator.get('nodedist',1)}

% if model.get('coord_system', 0) > 0:
m.cood_system     = ${model.get(['coord_system'])} 
% endif
