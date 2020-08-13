  exit_on_error=${model.get('exit_on_error', 'false')}
  exit_on_end=${model.get('exit_on_end')}
  verbosity=${model.get('verbosity', 2)}

model = '${model.get(['name'])}'
load_model(model)

