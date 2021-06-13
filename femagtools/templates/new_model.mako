
  exit_on_error=${model.get('exit_on_error', 'false')}
  exit_on_end=${model.get('exit_on_end', 'false')}
  verbosity=${model.get('verbosity', 1)}

model = '${model.get('name')}'
description = '${model.get('description','')}'
% if model.get('scratch_mode', 0):
new_model(model, '', 'scratch')
% else:
new_model_force(model, description)
%endif
