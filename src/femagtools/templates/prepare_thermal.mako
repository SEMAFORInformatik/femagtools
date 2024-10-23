-- prepare thermal model

save_model('cont')

function get_xy(elkey)
  -- get element xy
  local x, y
  x, y = get_elem_data ( "xycp", elkey)
  return x, y
end

function get_perm(elkey)
  -- get element permeability
  local x, y
  x, y = get_elem_data ( "perm", elkey)
  return x, y
end

function get_mag(elkey)
  -- get element rem
  local x, y
  x, y = get_elem_data ("mag", elkey)
  return x, y
end

function color_spel(spel_key, color)
  local el
  el = get_spel_data('elkeys', spel_key)
  for i = 1, #el do
     draw_elem(el[i], color, 0)
  end
end

get_spel_keys("spel_keys")      -- Get all superelements of the model

copper_spel = {}  -- red, 1
stator_lam_spel = {} -- darkblue, 9
rotor_lam_spel = {} -- darkblue, 9
shaft_spel = {} ---- lightgrey, 16
magnet_spel = {} -- black, 8
outer_air = {} -- darkred, 10
inner_air = {} -- lightgrey, 16
rotor_air = {} 
stator_air = {} -- black, 8
airgap = {} -- white, 8
magnet_pocket_spel = {}  -- red, 1
non_uniform_airgap = {} -- white, 8

-- search all components
for i =1,  #spel_keys do
  -- loop over all spel keys
  is_winding = get_spel_data("condtyp", spel_keys[i])
  els = get_spel_data('elkeys', spel_keys[i])
  mag1, mag2 = get_mag(els[1])
  perm1, perm2 = get_perm(els[1])
  xc, yc = get_xy(els[1])
  rc, pc = c2pd(xc, yc)

  if (is_winding == 0) then
     mctype = get_spel_data("mcvtyp", spel_keys[i])
     if (mctype == 0) then -- no curve
        -- check if it is air or magnet
        if math.abs(mag1) > 0 or math.abs(mag2) > 0 then
           -- is magnet
           table.insert(magnet_spel, spel_keys[i])
           def_mat_therm(xc, yc,'red', magn_density, magn_thcond, magn_thcap, 1)
           color_spel(spel_keys[i], 8) -- black magnet
        else
           -- is air
           if rc > dy1/2 then
              -- outer air
              table.insert(outer_air, spel_keys[i])
%if model.get('htc_outer', 0):
              def_heat_transfer(xc,yc,yellow,${model['htc_outer']}, 1.0)
%endif
              color_spel(spel_keys[i], 10) 

           elseif rc > dy2/2 and rc < da2/2 then
              -- rotor air
              -- check is magnet pocket/air pocket/non_uniform airgap
              table.insert(rotor_air, spel_keys[i])
             --def_mat_therm(xc, yc, 'yellow', 1.12,0.026,1007, 1)
              --color_spel(spel_keys[i], 3) -- yellow air

           elseif rc > da1/2 and rc < dy1/2 then
              -- stator air
              table.insert(stator_air, spel_keys[i])
              def_mat_therm(xc, yc, 'yellow', 1.19,0.15,1007, 1)
              color_spel(spel_keys[i], 8) 

           elseif rc > da2/2 and rc < da1/2 then
              table.insert(airgap, spel_keys[i])
              def_mat_therm(xc, yc, 'yellow', 1.19,0.063,1007, 1)
              color_spel(spel_keys[i], 7) -- white air

--           elseif rc > dy2/2 then
              -- airgap
--              table.insert(airgap, spel_keys[i])
--              color_spel(spel_keys[i], 3) -- yellow air
--              def_mat_therm(xc, yc, 'yellow', 1.19,0.15,1007, 1)

            elseif rc < dy2/2 then
              -- check if shaft or inner air
              if x0_shaft ~= nil and x0_shaft ~= 0.0 then
                  table.insert(shaft_spel, spel_keys[i])
                  def_mat_therm(xc, yc, 'lightgrey', shaft_density, shaft_thcond, shaft_thcap, 1)
              else
                  -- is inner air
                  table.insert(inner_air, spel_keys[i])
%if model.get('htc_inner', 0):
                  def_heat_transfer(xc,yc,yellow,${model['htc_inner']}, 1.0)
%endif
                  def_new_sreg(xc,yc, 'INAR', 'yellow')
                  color_spel(spel_keys[i], 16) -- light grey inner air
              end
           end
        end
     else
        -- check if it is stator / rotor
        if rc > m.fc_radius then
           color_spel(spel_keys[i], 9) -- magenta stator
           table.insert(stator_lam_spel, spel_keys[i])
           def_mat_therm(xc, yc, 'darkblue', stator_density, stator_thcond, stator_thcap, 1)
        else
           color_spel(spel_keys[i], 9) -- violet rotor
           table.insert(rotor_lam_spel, spel_keys[i])
           def_mat_therm(xc, yc, 'darkblue', rotor_density, rotor_thcond, rotor_thcap, 1)
        end
     end

  else
     table.insert(copper_spel, spel_keys[i])
     def_mat_therm(xc, yc, 'green', conductor_density, conductor_thcond, conductor_thcap, 1)
     color_spel(spel_keys[i], 1) --

  end

end

function is_exist(arr, val)
   local i
   if #arr == 0 then 
      return false 
   end 
   for i = 1, #arr do 
      if arr[i] == val then 
         return true
      end 
   end 
   return false
end 
-- identify magnet air pocket 

for i = 1, #magnet_spel do 
    bkeys = get_spel_data ( "bndkeys", magnet_spel[i])
    for j = 1, #bkeys do 
         elks = get_node_data ( "elkeys", bkeys[j] )
         for k = 1, #elks do 
            p1, p2 = get_elem_data("perm", elks[k])
            sek = get_elem_data("sekey", elks[k])
            if p1 == 1.0 then 
                  if (not is_exist(magnet_pocket_spel, sek)) and (not is_exist(airgap, sek)) then  
                     table.insert(magnet_pocket_spel, sek)
                     els = get_spel_data('elkeys', sek)
                     xc, yc = get_xy(els[1])
                     def_mat_therm(xc, yc, 'red', 1.12,0.026,1007, 1)
                     --table.remove(rotor_air, sek)
                     color_spel(sek, 1)
                  end 
            end 
         end 
    end 
end 



for i = 1, #rotor_air do 
  if rotor_air[i] ~= nil then 
    bkeys = get_spel_data ( "bndkeys", rotor_air[i])
    for j = 1, #bkeys do 
      x, y = get_node_data ( "xy", bkeys[j] )
      r, phi = c2pd(x, y) 
      if math.abs(r - da2/2) < 1e-5 then 
        if not is_exist(non_uniform_airgap,rotor_air[i]) then 
          table.insert(non_uniform_airgap, rotor_air[i])
          els = get_spel_data('elkeys', rotor_air[i])
          xc, yc = get_xy(els[1])
          def_mat_therm(xc, yc, 'yellow', 1.19,0.063,1007, 1)
          color_spel(rotor_air[i], 7)
        end 
      end 
    end 
  end 
end 

for i = 1, #rotor_air do 
  if rotor_air[i] ~= nil then 
    if is_exist(non_uniform_airgap, rotor_air[i]) then 
      rotor_air[i] = nil
    end    
  end 
end 

for i = 1, #rotor_air do 
   if is_exist(magnet_pocket_spel, rotor_air[i]) then 
     rotor_air[i] = nil
   end 
 end 
 

for i = 1, #rotor_air do 
   if rotor_air[i] ~= nil then 
      els = get_spel_data('elkeys', rotor_air[i])
      xc, yc = get_xy(els[1])
      def_mat_therm(xc, yc, 'red', 1.12,0.026,1007, 1)
   end 
end 

save_metafile(model..'.ps')
