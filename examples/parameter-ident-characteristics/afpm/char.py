"""
  Speed characteristics of Axialflux PM
  Ronald Tanner
"""
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
import femagtools.plot
import femagtools.machine
from femagtools.machine.utils import betai1


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    with open('dqpar.json') as fp:
        dqpars = json.load(fp)

    temp = [90, 90]
    pm = femagtools.machine.create_from_eecpars(temp, dqpars)

    i1max = 80
    Udc = 800
    u1max = round(0.9*Udc/np.sqrt(2)/np.sqrt(3))

    rep = []
    w1, T = pm.w1_imax_umax(i1max, u1max)
    iq, id = pm.iqd_tmech(T, w1/2/np.pi/pm.p)
    uq, ud = pm.uqd(w1, iq, id)
    u1 = np.linalg.norm((ud, uq))/np.sqrt(2)
    beta, i1 = betai1(iq, id)
    rep.append((w1/2/np.pi*60/pm.p, T, i1, beta/np.pi*180, u1))

    nmax = 4000/60
    w1 = 2*np.pi*nmax*pm.p
    iq, id, tq = pm.iqd_imax_umax(i1max, w1, u1max, T)
    beta, i1 = betai1(iq, id)
    uq, ud = pm.uqd(w1, iq, id)
    u1 = np.linalg.norm((ud, uq))/np.sqrt(2)
    rep.append((w1/2/np.pi*60/pm.p, tq, i1, beta/np.pi*180, u1))

    print("""
  n/rpm  T/Nm I1/A beta/°  U1/V
  -------------------------------
  """)
    for r in rep:
        print(f" {r[0]:5.0f}  {r[1]:5.1f} {r[2]:5.1f} {r[3]:5.1f} {r[4]:5.1f}")
    print()

    r = pm.characteristics(T, nmax, u1max)
    fig = femagtools.plot.characteristics(r)
    fig.savefig('speedchar.pdf')
