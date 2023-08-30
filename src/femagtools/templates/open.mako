  exit_on_error=${model.get('exit_on_error', 'false')}
  exit_on_end=${model.get('exit_on_end')}
  verbosity=${model.get('verbosity', 1)}

model = '${model.get(['name'])}'
load_model(model)

m.num_poles = get_dev_data("num_poles")
m.npols_gen = get_mod_data("num_poles")
m.arm_length = get_dev_data("arm_length")
