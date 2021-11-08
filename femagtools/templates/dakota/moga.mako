# Dakota Input Template File

environment
  tabular_data
    tabular_data_file = 'femag.dat'

method
  moga
% if set(('num_generations', 'percent_change')).intersection(study.keys()):
  convergence_type metric_tracker
% if 'num_generations' in study:  
    num_generations = ${study['num_generations']}
% endif
% if 'percent_change' in study:  
    percent_change = ${study['percent_change']}
% endif
% endif
  #final_solutions = 3
  print_each_pop
#    fitness_type domination_count
#    niching_type distance 0.05 0.05
#    replacement_type below_limit = 6
#    replacement_type roulette_wheel
#    postprocessor_type orthogonal_distance 0.05 0.05     
#    max_function_evaluations = ${study.get('max_function_evaluations', 1000)}

variables
% if [d['bounds'] for d in study['decision_vars']]:
  continuous_design = ${len(study['decision_vars'])}
    lower_bounds  ${' '.join([str(d['bounds'][0]) for d in study['decision_vars']])}
    upper_bounds  ${' '.join([str(d['bounds'][1]) for d in study['decision_vars']])}
    descriptors   ${' '.join([f"\"{d['name']}\"" for d in study['decision_vars']])}
% endif

interface
  id_interface = 'FEMAG'
  analysis_driver = 'python -m femagtools.dakota_femag'
    fork batch 
    parameters_file = 'params.in'
    results_file    = 'results.out'
    file_save
    
responses
  objective_functions = ${len(study['objective_vars'])}
  descriptors  ${' '.join([f"\"{o['name']}\"" for o in study['objective_vars']])}
  # sense 'max' 'min' ..
  no_gradients
  no_hessians
