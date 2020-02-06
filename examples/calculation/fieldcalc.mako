
-- demo template for custom calculations
---
maxit=300
du_u0=1e-3

m.num_par_wdgs = ${model.get('num_par_wdgs', 1)}
phi = ${model.get('phi')}
I1 = ${model.get('current')}*math.sqrt(2.0)/m.num_par_wdgs
      k = 0
      for ii=1,3 do
        cur = I1*math.cos(phi/180.0*math.pi+k*math.pi/3.0)
        def_curr_wdg(1,cur) 
        k = k+2
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
