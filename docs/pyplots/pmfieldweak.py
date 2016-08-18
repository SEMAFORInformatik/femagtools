import numpy as np
import matplotlib.pyplot as pl
import femagtools

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

tq = 170.0
u1 = 340.0

iqx, idx = pm.iqd_torque(tq)
w1 = pm.w1_u(u1, idx, iqx)
i1 = np.linalg.norm(np.array((iqx, idx)))

betaopt = np.arctan2(idx, iqx)/np.pi*180
print("N {0:8.1f} I1 {1:8.1f} Beta {2:4.1f}".format(w1/2/np.pi, i1, betaopt))
       
fig, ax = pl.subplots()

id = np.linspace(i1*np.sin(beta[0]/180*np.pi),
                 i1*np.sin(beta[-1]/180*np.pi))
iq = [pm.iq_u(w1, u1, ix) for ix in id]
ax.plot(id, iq, label='U1={} V'.format(u1))

for tq in (tq,): # (50, 100, 150, 200, 250):
    i1x = pm.i1_torque(tq, beta[0]/180*np.pi)
    iqmin = i1x*np.cos(beta[0]/180*np.pi)
    i1x = pm.i1_torque(tq, beta[-1]/180*np.pi)
    iqmax = i1x*np.cos(beta[-1]/180*np.pi)
    iq = np.linspace(iqmin, iqmax, 20)
    id = np.array([pm.id_torque(tq, ix) for ix in iq])
    ax.plot(id, iq, label='Tq={} Nm'.format(tq))
    #mini = np.linalg.norm(np.vstack((id, iq)), axis=0).argmin()
    #ax.plot(id[mini], iq[mini], 'ro')
    #iqx, idx = pm.iqd_torque(tq)
    #ax.plot(idx, iqx, 'bo')

#ax.plot(idx, iqx, 'ro')
ax.annotate('f1={0:3.0f} Hz'.format(w1/np.pi/2),
            xy=(idx, iqx), xytext=(1.1*idx, 1.5*iqx),
            arrowprops=dict(arrowstyle="->"))

ax.arrow(0, 0, idx+0.075*i1, iqx-0.08*i1, color='r',
         head_width=0.05*i1, head_length=0.08*i1)
ax.text(1.1*idx, 0.4*iqx, r'$I_1={0:4.1f} A$'.format(i1), fontsize=18)
ax.arrow(0, 0, 0, 170, color='k', head_width=0.05*i1, head_length=0.08*i1)

ax.annotate("",
            xy=(0.36*idx, 0.36*iqx), xycoords='data',
            xytext=(0, 0.4*i1), textcoords='data',
            arrowprops=dict(arrowstyle="->", #linestyle="dashed",
                            color="0.5",
                            shrinkA=5, shrinkB=5,
                            patchA=None,
                            patchB=None,
                            connectionstyle="arc3,rad=0.3"))
ax.text(0.52*idx, 0.58*iqx, r'$\beta={0:3.1f}^o$'.format(betaopt), fontsize=14)

# current limit
iqx0 = i1*np.cos(beta[0]/180*np.pi)
idx0 = i1*np.sin(beta[0]/180*np.pi)
id = np.linspace(idx0, idx)
#ax.plot(id, np.sqrt(i1**2 - id**2), label='I1 {0:8.1f} A'.format(i1))

# max speed
w1 = pm.w1_u(u1, idx0, iqx0)
iq = [pm.iq_u(w1, u1, ix) for ix in id]
#ax.plot(id, iq, label='n {0:8.1f} 1/min'.format(w1/2/np.pi/p*60))

print("max speed {}".format(w1/np.pi/2/p*60))
tq = pm.torque_iqd(iqx0, idx0)
i1x = pm.i1_torque(tq, beta[0]/180*np.pi)
iqmin = i1x*np.cos(beta[0]/180*np.pi)
i1x = pm.i1_torque(tq, beta[-1]/180*np.pi)
iqmax = i1x*np.cos(beta[-1]/180*np.pi)
iq = np.linspace(iqmin, iqmax, 20)
id = np.array([pm.id_torque(tq, ix) for ix in iq])
#ax.plot(id, iq, label='Tq {0:8.1f} Nm'.format(tq))
print("iq tab")

ax.grid()
ax.set_xlabel('Id / A')
ylim = ax.get_ylim()
ax.set_ylim([0, ylim[1]])
ax.set_ylabel('Iq / A')
legend = ax.legend(loc='upper center')

ax.set_aspect('equal')

pl.show()
