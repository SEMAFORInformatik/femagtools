
  for i=1,3 do
    def_curr_wdg(i,0.0)  
  end

  calc_field_single()

  post_models("induct(x)","b")    -- Calculate field distribution

  data=io.open("bag.dat","w")              -- Output in data file
--  data:write("phi [°]  Br [T]  Bph [T]\n");
  N = table.getn(b)                             -- Number of elements in array
  i = 1
  repeat
    data:write(string.format("%g %g %g\n",b[i],b[i+1],b[i+2]))
    i = i+3
  until i>=N
  io.close(data)                  -- Don't forget to close the file

 field_lines('field.svg',20)
 color_gradation( 0,0, tot, Babs, 0, 2.4, 'babs.svg')
