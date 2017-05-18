import numpy as np
import matplotlib.pyplot as pl
import femagtools.machine


#filename = '/home/tar/Documents/ScientificPython/download/TEST_001.BCH-1'
filename = 'TEST_001.BCH'
bch = femagtools.read_bchfile(filename)

beta = bch.ldq['beta']
i1 = bch.ldq['i1']

r10 = 0.1
ls = 1e-4
p = bch.machine['p']

#pm = femagtools.machine.PmRelMachineLdq(3, p,
#                                        r1=r10,
#                                        psid=bch.ldq['psid'],
#                                        psiq=bch.ldq['psiq'],
#                                        beta=beta,
#                                        i1=i1)

pm = femagtools.machine.PmRelMachineLdq(3, p,
                                        r1=r10, ls=ls,
                                        ld=bch.ldq['ld'],
                                        psim=bch.ldq['psim'],
                                        lq=bch.ldq['lq'],
                                        beta=beta,
                                        i1=i1)


Tmax = 170.0
u1 = 230.0
nmax = 70

iq, id = pm.iqd_torque(Tmax)
w1 = pm.w1_u(u1, iq, id)
w2 = pm.w2_i_u(femagtools.machine.betai1(iq, id)[1], u1)
#pmax = w1/p*Tmax
print("n1 = {} n2 = {}".format(60*w1/2/np.pi/p, 60*w2/2/np.pi/p))

r = pm.characteristics(Tmax, nmax, u1)
n = np.array(r['n'])
fig, axs = pl.subplots(2, sharex=True)
#fig.suptitle('PM-Machine Characteristics')

axs[0].plot(60*n, r['T'], 'b-', label='Torque')
axs[0].set_ylabel("Torque / Nm")
#ax[0].set_xticks(ticks)
axs[0].grid()
axs[0].legend(loc='upper right')
ax2 = axs[0].twinx()
ax2.plot(60*n, r['i1'], 'r-', label='Current')
ax2.set_ylabel("Current / A")
ax2.legend(loc='center right')

axs[1].plot(60*n, r['u1'], 'b-', label='Voltage')
axs[1].set_ylabel("Voltage / V",)
axs[1].set_xlabel("Speed / rpm")
axs[1].grid()
axs[1].legend(loc='upper left')
ax3 = axs[1].twinx()
ax3.plot(60*n, r['cosphi'], 'r-', label='Cos Phi')
ax3.set_ylabel("Cos Phi")
ax3.legend(loc='center right')
#pl.legend(('T', 'i1', 'u1'), loc='lower right')
#pl.ylim(0, 1.2)
#pl.xlabel('Speed [rpm]')
#pl.grid()

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
