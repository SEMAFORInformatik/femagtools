-- elkey: element key
-- magn:  magnetization (remanence [T])
-- alphaM: Angle of magnetization
-- murpm: Relative permeability
-- Bd/Hd: flux density / field strength direct to magnetization [T], [kA/m]
-- Bq/Hq: flux deensity field strength perpendicular to magnetization
-- alpha_h: Angle between magnetization and field strength vector

-- FEMAG build in
function dmg(elkey, magn_a, alpha_m, murpm, Bd, Bq, Hd, Hq, alpha_h, Hc_lim)
  local alpha_mn = alpha_m
  local magn_an = magn_a
  local mupm = mu0 * murpm
  if Hd < Hc_lim then
    magn_an = Bd - Hc_lim * 1e3 * mupm -- approx. the remanence
    if magn_an < 1.1e-5 then           -- avoid reversal of magn.
      magn_an = 1.1e-5
    end
  end
  iret = 1
  return magn_an, alpha_mn, iret
end
