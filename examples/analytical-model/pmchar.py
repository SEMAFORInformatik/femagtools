import numpy as np
import matplotlib.pyplot as pl
import femagtools.machine
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

filename = '../bch-erg/TEST_001.BCH'
bch = femagtools.read_bchfile(filename)

beta = bch.ldq['beta']
i1 = bch.ldq['i1']

r10 = 0.1
ls = 1e-4
p = bch.machine['p']

pm = femagtools.machine.PmRelMachineLdq(3, p,
                                        r1=r10, ls=ls,
                                        psid=bch.ldq['psid'],
                                        psiq=bch.ldq['psiq'],
                                        beta=beta,
                                        i1=i1)

#pm = femagtools.machine.PmRelMachineLdq(3, p,
#                                        r1=r10, ls=ls,
#                                        ld=bch.ldq['ld'],
#                                        psim=bch.ldq['psim'],
#                                        lq=bch.ldq['lq'],
#                                        beta=beta,
#                                        i1=i1)

Tmax = 170.0
u1 = 230.0
nmax = 4200/60

r = pm.characteristics(Tmax, nmax, u1)

n = np.array(r['n'])
fig, axs = pl.subplots(2, sharex=True)

axs[0].plot(60*n, r['T'], 'b-', label='Torque')
axs[0].set_ylabel("Torque / Nm")
axs[0].grid()
axs[0].legend(loc='upper right')
ax2 = axs[0].twinx()
ax2.plot(60*n, r['i1'], 'r-', label='Current')
ax2.set_ylabel("Current / A")
ax2.legend(loc='lower left')

axs[1].plot(60*n, r['u1'], 'b-', label='Voltage')
axs[1].set_ylabel("Voltage / V",)
axs[1].set_xlabel("Speed / rpm")
axs[1].grid()
axs[1].legend(loc='upper left')
ax3 = axs[1].twinx()
ax3.plot(60*n, r['cosphi'], 'r-', label='Cos Phi')
ax3.set_ylabel("Cos Phi")
ax3.legend(loc='lower right')

fig.tight_layout()
pl.show()

#pl.plot(60*n,np.array(T)/Tmax)
#pl.plot(60*n,np.array(r['i1'])/max(r['i1']))
#pl.plot(60*n,np.array(r['u1'])/max(r['u1']))
#pl.legend( ('T', 'i1', 'u1'), loc='lower right' )
#pl.ylim(0,1.2)
#pl.xlabel('Speed [rpm]')
#pl.grid()
#pl.show()
