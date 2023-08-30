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
Abar = get_sreg_data("area",rotorbar)

maxit = m.num_nonl_it
permode = 'restore'
du_u0 = m.error_perm/5

psi={}
L={}

% if model.get('bar_len', 0):
barlen = ${model.get('bar_len')*1e3}
% else:
p = m.num_poles/2
Dr = da2-m.slot_height
barlen = m.arm_length+math.pi*Dr/Q2/math.sin(math.pi*p/Q2)
% endif

state_of_problem('mag_dynamic')
file_bar = io.open("bar.dat","w")
dfreq = 8
f0 = 2
bag = {}
pos = {}
leakind = {}
currdens={1, 2, 3, 4} -- factor for current variation (Current densities A/mm²)
Nfreqs = 15 --
for i= 1, Nfreqs do
 freq = f0 + (i-1)*dfreq
 for k = 1, #currdens do
    cur = Abar*currdens[k]
    def_curr_wdg(rotorbar, cur, 0)

    calc_field_single({maxit=maxit, maxcop=du_u0,
        permode=permode, freq=freq})
      if k==1 then
        u_re, u_im = get_wdg_data("volt", rotorbar)
        rbar = u_re*barlen/cur
      end
      i_re, i_im = get_wdg_data("cur", rotorbar)
      flx1_re, flx1_im = get_wdg_data("flux", rotorbar)
      flx2_re, flx2_im = get_wdg_data("flux", stator)

      print(string.format("%g: %g %g %g",
         freq, i_re, flx1_re*m.arm_length, flx2_re*m.arm_length))

      leakind[k] = (flx1_re-flx2_re)*m.arm_length/cur
  end
  file_bar:write(string.format("%g %g ",
      freq, rbar))
  for k=1, #currdens do
    file_bar:write(string.format("%g ", leakind[k]))
  end
  file_bar:write("\n")
end
file_bar:close()
