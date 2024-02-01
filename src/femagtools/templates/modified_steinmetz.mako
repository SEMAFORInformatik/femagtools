
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
  alpha=${model['alpha']}              -- Hysteresis Frequency Coefficient Alpha
  beta=${model['beta']}              -- Hysteresis Frequency Coefficient Beta
  spweight= ${model['special_weight']}   -- Specific Weight Iron [gr/m3]
  fillfact=${model['fillfac']}           -- Fillfactor Iron <= 1

  -- Modified Steinmetz Iron Loss Model
  -- pvfe = ch*hi*(b**(alpha + beta*b)) + cw*(hi**2)*(b**2)
  hxx = Bx/fillfact                 -- Bx
  hyy = By/fillfact                 -- By
  b21 = math.sqrt(hxx*hxx+hyy*hyy)

  b = b21/basind
  hi = fnu/basfrq
  coeff = alpha+beta*b
  ph = kh*spweight*km*ch*hi*(b^coeff)   -- [W/m3]
  pw = spweight*km*cw*(b^2)*(hi^2)   -- [W/m3]
  pe = 0.0

  iret=1
  return ph, pw, pe, iret
end