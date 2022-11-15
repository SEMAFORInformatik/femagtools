-- eddy current simulation of a rotor bar
--
-- set rotor bar
xb, yb = 0, da2/2 - m.slot_height/2 - m.slot_h1
muer = 1
rel = 100
cur = {1, 0}
sigma1, sigma2 = get_dev_data( "cond_conduct" )
rotorbar = def_new_bar(xb,yb, "green", "U", cur[1],cur[2], "wi",
                       sigma2, muer, rel, "polar", 0, 0)
-- A = get_sreg_data("area",rotorbar)
-- print(A)

maxit = m.num_nonl_it
permode = 'restore'
du_u0 = m.error_perm/5

psi={}
L={}

Dr = da2-m.slot_height
% if model.get('bar_len', 0):
m.num_phases = 3 -- TODO: fix this
barlen = ${model.get('bar_len')*1e3}
% else:
p = m.num_poles/2
barlen = m.arm_length+math.pi*Dr/Q2/math.sin(math.pi*p/Q2)
% endif
--lfe = 1.4*m.arm_length -- get_dev_data("arm_length")

state_of_problem('mag_dynamic')
file_bar = io.open("bar.dat","w")
dfreq = 8
f0 = 2
Nfreqs = 15 --
for i= 1, Nfreqs do
 freq = f0 + (i-1)*dfreq
 calc_field_single({maxit=maxit, maxcop=du_u0,
        permode=permode, freq=freq})
  u_re, u_im = get_wdg_data("volt", rotorbar)
  i_re, i_im = get_wdg_data("cur", rotorbar)
  print(string.format("%g: %g %g", freq, u_re*barlen, u_im*barlen))
    file_bar:write(string.format("%g %g %g\n",
      freq, u_re*barlen, u_im*barlen))
end
file_bar:close()
