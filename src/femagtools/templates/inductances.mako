% if model.get('calc_inductances',0):

L={}
curr = 1/a
ksym = m.num_poles/m.npols_gen
curr = {
  {1/a, 0, 0},
  {0, 1/a, 0},
  {0, 0, 1/a} }

data=io.open("inductances.dat","w+")
for i=1,m.num_phases do
  for k=1,m.num_phases do
    def_curr_wdg(k,curr[i][k],0)
  end
  calc_field_single({
        maxit=1,permode='actual',freq=m.frequency})

  for k=1,m.num_phases do
    psi_Re,psi_Im=flux_winding_wk(k)
    data:write(string.format("%g ",
               math.sqrt(psi_Re^2+psi_Im^2)*m.arm_length/a*ksym))
  end
  data:write('\n')
end
io.close(data)
% endif
