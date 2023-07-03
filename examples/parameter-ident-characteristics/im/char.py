"""
  speed characteristics of IM
  Ronald Tanner
"""
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
import femagtools.plot
import femagtools.machine

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')


with open('impar.json') as fp:
    impars = json.load(fp)

temp = [90, 80]
im = femagtools.machine.create_from_eecpars(temp, impars)

T = 57
u1max = 230
nmax = 100

print("""
  n/rpm   T/Nm   I1/A  U1/V  cosphi  f/Hz
-----------------------------------------
""")

wm = im.wmfweak(u1max, im.psiref, T)
for tq, n in zip([T, 32.9, 20],
                 [wm/2/np.pi, 1450/60, 2900/60]):
    r = im.operating_point(T, n, u1max)
    print(f" {n*60:6.0f}  {tq:5.1f}  {r['i1']:5.1f}  {r['u1']:4.1f} {r['cosphi']:6.3f}  {r['f1']:5.1f}")
print()

r = im.characteristics(T, nmax, u1max)
fig = femagtools.plot.characteristics(r)
fig.savefig('speedchar.pdf')
