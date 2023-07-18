
-- no load single calculation 

maxit=300
du_u0=1e-3

get_wdg_keys("wkeys")

-- check the machine type PSM/ESM
magn = get_dev_data("magn_remanence")

if magn > 0 then -- IPM
  -- set cur to zero 
  for i = 1, #wkeys do   
    def_curr_wdg(wkeys[i],0, 0) 
  end 
else -- ESM
  for i = 1, #wkeys - 1 do   
    def_curr_wdg(wkeys[i],0, 0) 
  end 
end 
   
calc_field_single(maxit, reset, du_u0)

post_models("induct(x)","b")    -- Calculate field distribution

  data=io.open("bag.dat","w")              -- Output in data file
  N = table.getn(b)                             -- Number of elements in array
  i = 1
  repeat
    data:write(string.format("%g %g %g\n",b[i],b[i+1],b[i+2]))
    i = i+3
  until i>=N
  io.close(data)

color_gradation(0,0,"tot","Babs",0,0,"")
-- experimental (new femag-classic needed)
-- without grf_clear, overlay fieldlines with color gradation
add_field_lines("field.svg", 25)