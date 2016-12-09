% for p in model.get('plots', []):
  % if p == 'field_lines':
  field_lines('field.svg', 20)
  % else:
  color_gradation( 0,0, 'tot', '${p}', 0, 0, '${p}'..'.svg')
  % endif
% endfor
