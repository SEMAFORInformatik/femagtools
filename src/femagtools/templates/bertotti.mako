
function pfe(Bx,By,kh,fnu,km,z,d)

  -- Bx   flux density x component
  -- By   flux density y component
  -- kh   Hysterese Factor 
  -- fnu  frequency
  -- km   material factor
  -- ph   return for hysterese losses
  -- pw   return for eddy current losses
  -- iret status 

  -- Parameter

  basfrq=${model['base_frequency']}      -- Base Frequency for ch and cw [Hz]
  basind=${model['base_induction']}      -- Base Induction (Peak) [T]
  ch=${model['ch']}                      -- Fe Hysteresis Coefficient ch [W/kg]
  cw=${model['cw']}                      -- Fe Eddy Current Coefficient cw [W/kg]
  ce=${model['ce']}                      -- Fe Excess Coefficient cw [W/kg]
  hyscoef=${model['alpha']}              -- Hysteresis Frequency Coefficient
  spweight= ${model['special_weight']}   -- Specific Weight Iron [gr/m3]
  fillfact=${model['fillfac']}           -- Fillfactor Iron <= 1

  -- Bertotti Iron Loss Model
  -- pvfe = ch*f*(b**alpha) + cw*f**2*(b**2) + ce*f**1.5*(b**1.5)
  hxx = Bx/fillfact                 -- Bx
  hyy = By/fillfact                 -- By
  b21 = math.sqrt(hxx*hxx+hyy*hyy)
  b = b21/basind

  hi = fnu/basfrq
  hcw = hi^2
  hce = hi^1.5

  ph = kh*spweight*km*ch*hi*(b^hyscoef)   -- [W/m3]
  pw = spweight*km*cw*hcw*(b^2)     -- [W/m3]
  pe = spweight*km*ce*hce*(b^1.5)
  iret=1
  return ph, pw, pe, iret
end