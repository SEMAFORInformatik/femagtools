
-- no load single calculation 

maxit=300
du_u0=1e-3

get_wdg_keys("wkeys")

-- set cur to zero 
for i = 1, #wkeys do   
  def_curr_wdg(wkeys[i],0, 0) 
end 
% if model.get('noload_ex_cur', 0):
def_curr_wdg(wkeys[#wkeys], ${model['noload_ex_cur']}, 0) 
% endif

calc_field_single(maxit, reset, du_u0)

post_models("induct(x)","b_airgap")    -- Calculate field distribution

  data=io.open("bag.dat","w")              -- Output in data file
  N = table.getn(b_airgap)                             -- Number of elements in array
  i = 1
  repeat
    data:write(string.format("%g %g %g\n",b_airgap[i],b_airgap[i+1],b_airgap[i+2]))
    i = i+3
  until i>=N
  io.close(data)

-- experimental (new femag-classic needed)
export_calc_results('fieldcalc.vtu')