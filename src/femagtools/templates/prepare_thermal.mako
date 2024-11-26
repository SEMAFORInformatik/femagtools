-- prepare thermal model

save_model('cont')
state_of_problem("therm_static")
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
slot_opening = {}
slot_sreg_exist = 0
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
              key_exist_sl = get_sreg_key(xc, yc)
              if key_exist_sl <= 0 and slot_sreg_exist == 0 then
                  def_new_sreg(xc, yc, 'Slot', 'yellow')
                  slot_sreg_exist = 1
              else
                  add_to_sreg(xc, yc, 'Slot')
              end
           elseif rc > da2/2 and rc < da1/2 then
              table.insert(airgap, spel_keys[i])
              def_mat_therm(xc, yc, 'yellow', 1.19,0.15,1007, 1)
              add_to_sreg(xc, yc, 'Slot')
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
          def_mat_therm(xc, yc, 'yellow', 1.19,1.15,1007, 1)
          add_to_sreg(xc, yc, 'Slot')
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



---------------------------------------------
-- material slot insulation
---------------------------------------------

function get_boundary_node(num_slots)

   local ctr
   -- get the number of winding
   get_wdg_keys("wkeys")
   -- search subregions
   local srkey_sl = {}
   ctr = 1
   for i = 1, #wkeys do
      tmpksl = get_wdg_data("srkeys", wkeys[i])
      for j = 1, #tmpksl do
         srkey_sl[ctr] = tmpksl[j]
         ctr = ctr + 1
      end
   end
   -- search superelements
   local sekey_sl = {}
   ctr = 1
   for k = 1, #srkey_sl do
      sl_se = get_sreg_data("sekeys", srkey_sl[k])
      for kk =1, #sl_se do
         sl_nd = get_spel_data("bndkeys", sl_se[kk])
         x, y = get_node_data("xy", sl_nd[1])
         r, phi = c2pd(x, y)
         if phi < 360/num_slots then
            sekey_sl[ctr] = sl_se[kk]
            ctr = ctr + 1
         end
      end
   end

   local xn  = {}
   local yn = {}
   local node = {}
   local bndnodes = {}
   local bnd_unique = {}

   ctr = 1

   for i = 1, #sekey_sl do
      bnd = get_spel_data("bndkeys", sekey_sl[i])

      for j = 1, #bnd do
         bndnodes[ctr] = bnd[j]
         ctr = ctr + 1
         if bnd_unique[bnd[j]] == true then
            bnd_unique[bnd[j]] = false
         else
            bnd_unique[bnd[j]] = true
         end
      end
   end

   ctr =  1
   for j = 1, #bndnodes do
     x, y = get_node_data("xy", bndnodes[j])

     r, phi = c2pd(x, y)
     if (phi < 360/num_slots/2 - 0.05) and (bnd_unique[bndnodes[j]] == true)  then

        node[ctr] = bndnodes[j]
        xn[ctr] = x
        yn[ctr] = y
        ctr = ctr + 1
       end
   end


   local indx = {1, math.floor(#node/4), math.floor(#node/2), math.floor(#node/4*3), #node}

   local x_new = {}
   local y_new = {}

   for i = 1, 10 do
      x_new[i] = 0
      y_new[i] = 0
   end

   for i = 1, #indx do
      x_new[i] = xn[indx[i]]
      y_new[i] = yn[indx[i]]
      r, phi = c2pd(x_new[i], y_new[i])
      x_new[10 - (i-1)], y_new[10 - (i-1)] = pd2c(r, 360/num_slots - phi)
   end

   local rn = {}
   local phin = {}

   for i = 1, #x_new do
      rn[i], phin[i] = c2pd(x_new[i], y_new[i])
   end
   return rn, phin
end


rn, phin = get_boundary_node(m.tot_num_slot)
thickness = ${model.windings.get('slot_insulation_thickness', 0.15e-3)*1e3}
if thickness == 0.0 then 
   thickness = 0.15
end 
conductivity = ${model.get('slot_insul_cond', 0.31)}

for i = 1,m.num_sl_gen do
x = {}
y = {}

for j = 1, #rn do
   x[j], y[j] = pd2c(rn[j], phin[j] + (i-1)*360/m.tot_num_slot)
   point(x[j], y[j],"black","x")
end

def_insulation_by_nodechain(thickness,conductivity,
      x[1],y[1],
      x[2],y[2],
      x[3],y[3],
      x[4],y[4],
      x[5],y[5],
      x[6],y[6],
      x[7],y[7],
      x[8],y[8],
      x[9],y[9],
      x[10],y[10],
      x[1], y[1]
      )

end
save_metafile('insulation.ps')

state_of_problem("mag_static")
