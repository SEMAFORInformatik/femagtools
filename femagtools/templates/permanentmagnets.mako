

-- pm magnets
if tmp.mag_exists > 0 then
  alfa = tmp.mag_alpha
  for i=0, m.npols_gen-1 do
    for n=1, tmp.mag_exists do
      r, p = c2pr(tmp.xmag[n], tmp.ymag[n])
      phi = i*alfa+p
      x0, y0 = pr2c(r, phi)
      phi_orient = i*alfa+tmp.mag_orient[n]
      orient = phi_orient*180/math.pi
      if ( i % 2 == 0 ) then
        def_mat_pm(x0, y0, red, m.remanenc, m.relperm,
                   orient-180, m.parallel, 100)
      else
        def_mat_pm(x0, y0, green, m.remanenc, m.relperm,
                   orient, m.parallel, 100)
      end
      if tmp.mag_mirrored then
        phi = (i+1)*alfa-p
        x0, y0 = pr2c(r, phi)
        phi_orient = (i+1)*alfa-tmp.mag_orient[n]
        orient = phi_orient*180/math.pi
        if ( i % 2 == 0 ) then
          def_mat_pm(x0, y0, red, m.remanenc, m.relperm,
                     orient-180, m.parallel, 100)
        else
          def_mat_pm(x0, y0, green, m.remanenc, m.relperm,
                     orient, m.parallel, 100)
        end
      end
    end
  end
end