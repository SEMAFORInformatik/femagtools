-- Plots
% for p in model.get('plots', []):
% if not isinstance(p, str):
% if p[0] == 'field_lines':
% if len(p) < 2:
field_lines('field.svg', 20)
% else:
field_lines('field.svg', ${p[1]})
% endif
% else:
% if len(p) < 3:
color_gradation( 0,0, 'tot', '${p[0]}', 0, 0, '${p[0]}'..'.svg')
% else:
color_gradation( 0,0, 'tot', '${p[0]}', ${p[1]}, ${p[2]}, '${p[0]}'..'.svg')
% endif
% endif
% else:
% if p == 'field_lines':
field_lines('field.svg', 20)
% else:
color_gradation( 0,0, 'tot', '${p}', 0, 0, '${p}'..'.svg')
% endif
% endif
% endfor
