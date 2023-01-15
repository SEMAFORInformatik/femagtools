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


mass_stator = ${model['stator_material'][0]}
emodulus_st = ${model['stator_material'][1]*1e9}
possion_st = ${model['stator_material'][2]}

-- subregion stator yoke 
x,y = pd2c(dy1/2-1, 360/m.num_sl_gen/4)
def_temp_mech(x,y,"darkblue",20)
def_mat_mech(x,y,"darkblue",mass_stator, emodulus_st, possion_st,0,0,0)
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
def_mat_mech(x1,y1,"darkblue",mass_stator, emodulus_st, possion_st,0,0,0)

-- define material for the slots

get_sreg_keys("sreg_keys")

sreg_winding = {} 

ctr = 1
for i =1, #sreg_keys do 
   wind_key = get_sreg_data ( "wbkey", sreg_keys[i] )
   if wind_key > 0 then  
      sreg_winding[ctr] = sreg_keys[i]
      ctr = ctr + 1
   end 
end 

ctr = 1 
el_keys_winding = {}
for j = 1, #sreg_winding do 
   tmp_sekeys = get_sreg_data ( "sekeys", sreg_winding[j] )
   for k = 1, #tmp_sekeys do 
        tmp_elkeys = get_spel_data ("elkeys", tmp_sekeys[k] )
        el_keys_winding[ctr] = tmp_elkeys[1]
        ctr = ctr + 1
   end 
end 

mass_slot = ${model['slot_material'][0]}
emodulus_slot = ${model['slot_material'][1]*1e9}
possion_slot = ${model['slot_material'][2]}

for kk = 1, #el_keys_winding do 

  x0, y0 = get_elem_data ( "xycp", el_keys_winding[kk])
  def_temp_mech(x0, y0,"lightgrey",20)
  def_mat_mech(x0, y0,"lightgrey",mass_slot, emodulus_slot, possion_slot,0,0,0)

end 

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