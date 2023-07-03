"""
  speed characteristics of EESM
  Ronald Tanner
"""
import logging
import numpy as np
import json
import femagtools.plot
import femagtools.machine
from femagtools.machine.utils import betai1
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

with open('eecpars.json') as fp:
    eecpars = json.load(fp)

temp = [90, 90]
eesm = femagtools.machine.create_from_eecpars(temp, eecpars)

i1max = 275
udc = 400
u1max = 0.9*udc/np.sqrt(3)/np.sqrt(2)
w1, Tmax = eesm.w1_imax_umax(i1max, u1max)

speed = w1/2/np.pi/eesm.p
iq, id, iex = eesm.iqd_tmech(Tmax, speed)
uq, ud = eesm.uqd(w1, iq, id, iex)
u1 = np.linalg.norm((ud, uq))/np.sqrt(2)

print("""
  n/rpm  T/Nm   I1/A   Iex/A  U1/V
-----------------------------------
""")
for n, T in zip([2000/60, speed, 6000/60, 10000/60],
                [120, Tmax, 100, 60]):
    w1 = 2*np.pi*n*eesm.p
    iq, id, iex, tqx = eesm.iqd_tmech_umax(T, w1, u1max)
    beta, i1 = betai1(iq, id)
    uq, ud = eesm.uqd(w1, iq, id, iex)
    u1 = np.linalg.norm((ud, uq))/np.sqrt(2)
    print(f" {60*n:6.0f}  {T:5.1f}  {i1:5.1f}  {iex:4.1f}  {u1:5.1f}")
print()
nmax = 12000/60
r = eesm.characteristics(Tmax, nmax, u1max)
fig = femagtools.plot.characteristics(r)
fig.savefig('speedchar.pdf')
