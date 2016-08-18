import numpy as np
import matplotlib.pyplot as pl
import femagtools


def torque(T, pmax, wm):
    """shaft torque as a function of rotor angular speed"""
    if wm <= pmax / T:
        return T
    return pmax / wm

p = 4
r1 = 0.0806
le = 0.0
ls = 0.0
wind_temp = 20.0
ld = [0.0014522728, 0.0014522728]
lq = [0.0032154, 0.0038278836]
psim = [0.11171972000000001, 0.11171972000000001]
i1 = [80.0]
beta = [0.0, -41.1]

pm = femagtools.PmRelMachineLdq(3, p,
                                psim,
                                ld,
                                lq,
                                r1,
                                beta,
                                i1)

pmax = 60e3
Tmax = 170.0
u1 = 340.0
n = np.linspace(0, 75, 20)
T = [torque(Tmax, pmax, 2*np.pi*nx) for nx in n]

r = pm.characteristics(T, n, u1)

fig, axs = pl.subplots(2, sharex=True)
#fig.suptitle('PM-Machine Characteristics')
axs[0].plot(60*n, T, 'b-', label='Torque')
axs[0].set_ylabel("Torque / Nm")
#ax[0].set_xticks(ticks)
axs[0].grid()
axs[0].legend(loc='upper center')
ax2 = axs[0].twinx()
ax2.plot(60*n, r['i1'], 'r-', label='Current')
ax2.set_ylabel("Current / A")
ax2.legend(loc='lower right')

axs[1].plot(60*n, r['u1'], 'b-', label='Voltage')
axs[1].set_ylabel("Voltage / V",)
axs[1].set_xlabel("Speed / rpm")
axs[1].grid()
axs[1].legend(loc='upper center')
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
