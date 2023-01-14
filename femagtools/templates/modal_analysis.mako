
-- delete all windings
del_all_coils()

-- TODO: define material for the slots

<% 
import sys 
platform = sys.platform
%>
platform = "${platform}"

if platform == "win32" then 
    point_to_me("wmefemagw64.exe")
else
    point_to_me("mefemag64")
end 

-- subregion stator yoke 
x,y = pd2c(dy1/2-1, 360/m.num_sl_gen/4)
def_temp_mech(x,y,"darkblue",20)
def_mat_mech(x,y,"darkblue",${model['mass_density']}, ${model['youngs_modulus']*1e9},0.3,0,0,0)
-- subregion stator teeth 
x,y = pd2c(da1/2+1, 0)
key = get_elem_key(x, y)
-- identify if the area is slot or tooth
px, py = get_elem_data("perm", x, y)

if px > 1 then 
    x1 = x
    y1 = y
else
    x1,y1 = pd2c(da1/2+1, 360/m.num_sl_gen/2 - 0.5)
end 

def_temp_mech(x1,y1,"darkblue",20)
def_mat_mech(x1,y1,"darkblue",${model['mass_density']}, ${model['youngs_modulus']*1e9},0.3,0,0,0)

ev = ${model.get('num_modes', 20)}
calc_mech_eigenvalue(ev)

% if model['export_figure']:
export_fig = true
% else: 
export_fig = false
%endif

ctr = 1
for k = 1, ev, 1 do
    frequency = get_eigenfrequency (k)
    if frequency > 20 then
        if export_fig then 
	        draw_eigenform ( k,model .. "_Mode_" .. string.format('%03d', ctr) .. ".ps" )
        end 
	    iret = export_eigenvectors(1,k, model .. "_Mode_" .. string.format('%03d', ctr) .. ".txt")  -- mass normalized
        ctr = ctr + 1
    end

end


save_model('close')