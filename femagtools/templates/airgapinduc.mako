  post_models("induct(x)","bagx")    -- Calculate field distribution
  data=io.open("bag.dat","w")     -- Output in data file
  for i = 1, #bagx, 3 do
    data:write(string.format("%g %g %g\n",bagx[i],bagx[i+1],bagx[i+2]))
  end
  io.close(data)
