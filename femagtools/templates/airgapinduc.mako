  post_models("induct(x)","b")    -- Calculate field distribution

  data=io.open("bag.dat","w")              -- Output in data file
  N = table.getn(b)                             -- Number of elements in array
  i = 1
  repeat
    data:write(string.format("%g %g %g\n",b[i],b[i+1],b[i+2]))
    i = i+3
  until i>=N
  io.close(data)                  -- Don't forget to close the file
