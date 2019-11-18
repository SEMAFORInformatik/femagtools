  post_models("induct(x)","b")    -- Calculate field distribution
  data=io.open("bag.dat","w")     -- Output in data file
  for i = 1, table.getn(b), 3 do
    data:write(string.format("%g %g %g\n",b[i],b[i+1],b[i+2]))
  end
  io.close(data)
