"""
  femagtools.pm
  ~~~~~~~~~~~~~

  PM/Rel synchronous machine (SPM, IPM, RM) electrical circuit model

  Copyright 2022: Semafor Informatik & Energie AG, Switzerland
"""
import logging
import numpy as np
import numpy.linalg as la
from .utils import iqd, betai1, skin_resistance, dqparident
import scipy.optimize as so
import scipy.interpolate as ip

logger = logging.getLogger(__name__)


def parident(workdir, engine, temp, machine,
             magnetizingCurves, magnetMat, condMat,
             **kwargs):
    return dqparident(workdir, engine, temp, machine,
                      magnetizingCurves, magnetMat, condMat,
                      **kwargs)


class PmRelMachine(object):
    """Abstract base class for PmRelMachines

    Args:
        m: number of winding phases
        p: number of pole pairs
        r1: stator winding resistance (in Ohm)
        ls: leakage inductance in H
    """

    def __init__(self, m, p, r1, ls, **kwargs):
        self.p = p
        self.m = m
        self.r1 = r1
        self.ls = ls
        self.io = (1, -1)
        self.fo = 50.0
        self.tcu1 = 20
        self.zeta1 = 0.2
        self.gam = 0.7
        self.kh = 2
        for k in kwargs.keys():
            setattr(self, k, kwargs[k])

        self.plexp = {'styoke_hyst': 1.0,
                      'stteeth_hyst': 1.0,
                      'styoke_eddy': 2.0,
                      'stteeth_eddy': 2.0,
                      'rotor_hyst': 1.0,
                      'rotor_eddy': 2.0}
        self._losses = {k: lambda x, y: 0 for k in (
            'styoke_hyst', 'stteeth_hyst',
            'styoke_eddy', 'stteeth_eddy',
            'rotor_hyst', 'rotor_eddy',
            'magnet')}

    def rstat(self, w):
        """stator resistance"""
        return skin_resistance(self.r1, w, self.tcu1, self.zeta1,
                               self.gam, self.kh)

    def torque_iqd(self, iq, id):
        "torque at q-d-current"
        psid, psiq = self.psi(iq, id)
        tq = self.m*self.p/2*(psid*iq - psiq*id)
        return tq

    def iqd_torque(self, torque):
        """return minimum d-q-current for torque"""
        if np.abs(torque) < 1e-2:
            return (0, 0)
        res = so.minimize(lambda iqd: la.norm(iqd), self.io, method='SLSQP',
                          constraints=({'type': 'eq',
                                        'fun': lambda iqd:
                                        self.torque_iqd(*iqd) - torque}))
        if not res.success:
            raise ValueError(f'Torque {torque} out of current range')
        return res.x

    def tmech_iqd(iq, id, n, kpfe, pfw):
        """return shaft torque with given d-q current, iron loss correction factor
          and friction windage losses"""
        f1 = self.p*n
        plfe1 = self.iqd_plfe1(iq, id, f1)
        plfe2 = self.iqd_plfe2(iq, id, f1)
        plfe = kpfe * (plfe1 + plfe2)
        pmag = self.iqd_plmag(iq, id, f1)
        return self.torque_iqd(iq, id) - (plfe + pmag + pfw)/(2*np.pi*n)

    def uqd(self, w1, iq, id):
        """return uq, ud of frequency w1 and d-q current"""
        psid, psiq = self.psi(iq, id)
        uqd = (self.r1*iq + w1*(self.ls*id + psid),
               self.r1*id - w1*(self.ls*iq + psiq))
        logger.debug('beta i1 %s u1 %f', betai1(iq, id), la.norm(uqd))
        return uqd

    def w1_umax(self, u, iq, id):
        """return frequency w1 at given voltage u and id, iq current

        Keyword arguments:
        u -- the maximum voltage (RMS)
        iq, id -- the d-q currents"""
        w10 = np.sqrt(2)*u/la.norm(self.psi(iq, id))
        return so.fsolve(lambda w1:
                         la.norm(self.uqd(w1, iq, id))-u*np.sqrt(2), w10)[0]

    def w1_u(self, u, iq, id):
        """return frequency w1 at given voltage u and id, iq current
        (obsolete, use w1_umax)"""
        return self.w1_umax(u, iq, id)

    def w1max(self, u, iq, id):
        """return max frequency w1 at given voltage u and d-q current
        (obsolete, use w1_umax)"""
        return self.w1_umax(u, iq, id)

    def w2_imax_umax(self, imax, umax, maxtorque=True):
        """return frequency at max current and max voltage"""
        w, info, ier, mesg = so.fsolve(lambda x: np.linalg.norm(
            self.uqd(x, *iqd(-np.pi/2, imax))) - umax*np.sqrt(2),
            np.sqrt(2)*umax/la.norm(self.psi(*self.io)),
            full_output=True)
        if ier == 1:
            return w[0]

        logger.warn("w2_imax_umax ier=%d imax %f", ier, imax)
        raise ValueError("w2_imax_umax {} imax {}".format(mesg, imax))

    def beta_u(self, w1, u, i1):
        "beta at given frequency, voltage and current"
        return so.fsolve(lambda b:
                         la.norm(self.uqd(w1, *(iqd(b, i1))))-u*np.sqrt(2),
                         np.arctan2(self.io[1], self.io[0]))[0]

    def iq_u(self, w1, u, id):
        "iq at given frequency, voltage and id current"
        iq0 = max(self.io[0]/4, id*np.tan(self.betarange[0]))
        return so.fsolve(lambda iq:
                         la.norm(
                             self.uqd(w1, iq, np.array([id])))-u*np.sqrt(2),
                         iq0)[0]

    def iqd_uqd(self, w1, uq, ud):
        "return iq, id current at given frequency, voltage"
        return so.fsolve(lambda iqd:
                         np.array((uq, ud)) - self.uqd(w1, *iqd),
                         (0, self.io[1]))

    def i1_torque(self, torque, beta):
        "return i1 current with given torque and beta"
        i1, info, ier, mesg = so.fsolve(
            lambda i1: self.torque_iqd(*iqd(beta, i1))-torque,
            self.io[0],
            full_output=True)
        if ier == 1:
            return i1
        raise ValueError("no solution found for torque {}, beta {}".format(
            torque, beta))

    def i1_voltage(self, w1, u1, beta):
        "return i1 current with given w1, u1 and beta"
        i1, info, ier, mesg = so.fsolve(
            lambda i1: la.norm(self.uqd(w1, *iqd(beta, i1)))-np.sqrt(2)*u1,
            la.norm(self.io),
            full_output=True)
        if ier == 1:
            return i1
        raise ValueError("{} for w1 {}, u1 {}, beta {}".format(
            mesg, w1, u1, beta))

    def id_torque(self, torque, iq):
        "return d current with given torque and d-current"
        id0 = min(self.io[1]/4, iq/np.tan(self.betarange[0]))
        return so.fsolve(lambda id: self.torque_iqd(np.array([iq]), id)-torque, id0)[0]

    def iqd_torque_umax(self, torque, w1, u1max):
        "return d-q current and torque at stator frequency and max voltage"
        iq, id = self.iqd_torque(torque)
        # check voltage
        if la.norm(self.uqd(w1, iq, id)) <= u1max*np.sqrt(2):
            return (iq, id, torque)
        # decrease psi (flux weakening mode), let i1 == i1max
        iqd, info, ier, mesg = so.fsolve(
            lambda iqd: (la.norm(self.uqd(w1, *iqd)) - u1max*np.sqrt(2),
                         self.torque_iqd(*iqd) - torque),
            (iq, id),
            full_output=True)
        if ier != 1:  # didn't converge
            return self.mtpv(w1, u1max, betai1(iq, id)[1],
                             maxtorque=torque > 0)
        return iqd[0], iqd[1], self.torque_iqd(iq, id)

    def iqd_torque_imax_umax(self, torque, n, umax):
        """return iq, id, torque for constant torque or field weakening"""
        iq, id = self.iqd_torque(torque)
        w1 = 2*np.pi*n*self.p
        # Constant torque range
        if np.linalg.norm(self.uqd(w1, iq, id)) <= umax*np.sqrt(2):
            return (iq, id, torque)
        # Field weaking range
        imax = betai1(iq, id)[1]
        iq, id = self.iqd_imax_umax(imax, w1, umax, maxtorque=torque > 0)
        return iq, id, self.torque_iqd(iq, id)

    def iqd_imax_umax(self, i1max, w1, u1max, maxtorque=True):
        """return d-q current at stator frequency and max voltage
        and max current (for motor operation if maxtorque else generator operation)"""

        if maxtorque:
            btab = np.linspace(max(-np.pi/2, self.betarange[0]), 0)
        else:
            btab = np.linspace(max(-np.pi, self.betarange[0]), -np.pi/2)
        u1b = [np.linalg.norm(self.uqd(w1, *iqd(b, i1max)))
               for b in btab]
        if np.max(u1b) > np.sqrt(2)*u1max:
            beta0=btab[0]
            for bx, ux in zip(btab, u1b):
                if ux/np.sqrt(2)>u1max:
                    break
                beta0 = bx
            beta, info, ier, mesg = so.fsolve(
                lambda b: la.norm(
                    self.uqd(w1, *iqd(b, i1max))) - u1max*np.sqrt(2),
                beta0,
                full_output=True)
            if ier == 1:
                return iqd(beta[0], i1max)
        return self.mtpv(w1, u1max, i1max, maxtorque)[:2]

    def mtpa(self, i1):
        """return iq, id, torque at maximum torque of current i1"""
        sign = -1 if i1 > 0 else 1
        b0 = 0 if i1 > 0 else -np.pi
        bopt, fopt, iter, funcalls, warnflag = so.fmin(
            lambda x: sign*self.torque_iqd(*iqd(x, abs(i1))), b0,
            full_output=True,
            disp=0)
        iq, id = iqd(bopt[0], abs(i1))
        return [iq, id, sign*fopt]

    def mtpv(self, w1, u1, i1max, maxtorque=True):
        """return d-q-current, torque for voltage and frequency
        with maximum (maxtorque=True) or minimum torque """
        sign = -1 if maxtorque else 1
        i0 = (-sign*self.i1range[1]/10, -self.i1range[1]/10)
        res = so.minimize(
            lambda iqd: sign*self.torque_iqd(*iqd),
            i0, method='SLSQP',
            constraints=(
                {'type': 'ineq',
                 'fun': lambda iqd:
                 np.sqrt(2)*u1 - la.norm(self.uqd(w1, *iqd))},
                {'type': 'ineq',
                 'fun': lambda iqd:
                 i1max - betai1(*iqd)[1]}))

        return res.x[0], res.x[1], sign*res.fun

    def _set_losspar(self, pfe):
        self.fo = pfe['speed']*self.p
        ef = pfe.get('ef', [2.0, 2.0])
        hf = pfe.get('hf', [1.0, 1.0])
        self.plexp = {'styoke_hyst': hf[0],
                      'stteeth_hyst': hf[0],
                      'styoke_eddy': ef[0],
                      'stteeth_eddy': ef[0],
                      'rotor_hyst': hf[1],
                      'rotor_eddy': ef[1]}
        #                          'magnet'):

    def betai1_plcu(self, i1, w1=0):
        return self.m*self.rstat(w1)*i1**2

    def iqd_plcu(self, iq, id, w1=0):
        return self.m*self.rstat(w1)*(iq**2+id**2)/2

    def iqd_plcu1(self, iq, id, w1):
        return self.iqd_plcu(iq, id, w1)

    def iqd_plcu2(self, iq, id):
        return np.zeros(np.asarray(iq).shape)

    def betai1_losses(self, beta, i1, f):
        return np.sum([self.betai1_plfe1(beta, i1, f),
                       self.betai1_plfe2(beta, i1, f),
                       self.betai1_plmag(beta, i1, f),
                       self.betai1_plcu(i1), 2*np.pi*f], axis=0)

    def iqd_losses(self, iq, id, f):
        return np.sum([self.iqd_plfe1(iq, id, f),
                       self.iqd_plfe2(iq, id, f),
                       self.iqd_plmag(iq, id, f),
                       self.iqd_plcu(iq, id, 2*np.pi*f)], axis=0)

    def characteristics(self, T, n, u1max, nsamples=30):
        """calculate torque speed characteristics.
        return dict with list values of
        id, iq, n, T, ud, uq, u1, i1,
        beta, gamma, phi, cosphi, pmech, n_type

        Keyword arguments:
        T -- the maximum torque or the list of torque values in Nm
        n -- the maximum speed or the list of speed values in 1/s
        u1max -- the maximum voltage in V rms
        nsamples -- (optional) number of speed samples
        """
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        if np.isscalar(T):
            iq, id = self.iqd_torque(T)
            i1max = betai1(iq, id)[1]
            w1 = self.w1_umax(u1max, iq, id)
            w1max = self.w1_umax(u1max, *self.iqdmin(i1max))
            nmax = max(w1, w1max)/2/np.pi/self.p

            n1 = min(w1/2/np.pi/self.p, nmax)
            r['n_type'] = n1
            logger.info("Type speed %f n: %f nmax %f",
                        60*n1, 60*n, 60*nmax)
            try:
                w1 = self.w2_imax_umax(i1max, u1max, maxtorque=(T > 0))
                n2 = w1/2/np.pi/self.p
                iqmtpv, idmtpv, tq = self.mtpv(
                    w1, u1max, i1max, maxtorque=T > 0)
                if not self._inrange((iqmtpv, idmtpv)):
                    if n > 0:
                        n2 = min(nmax, n)
                    else:
                        n2 = nmax
                logger.info("n1: %f n2: %f ",
                            60*n1, 60*n2)
            except ValueError:
                if n > 0:
                    n2 = min(nmax, n)
                else:
                    n2 = nmax
            if n > 0:
                nmax = min(nmax, n)
                speedrange = sorted(
                    list(set([nx for nx in [n1, n2, n] if nx < 1.01*nmax])))
            else:
                speedrange = sorted(list(set([n1, n2])))
            logger.info("speedrange %s", speedrange)
            speedrange.insert(0, 0)
            n3 = speedrange[-1]
            nstab = [int(nsamples*(x1-x2)/n3)
                     for x1, x2 in zip(speedrange[1:],
                                       speedrange)]
            for nx in np.linspace(0, n1, nstab[0]):
                r['id'].append(id)
                r['iq'].append(iq)
                r['n'].append(nx)
                r['T'].append(T)

            n1 = speedrange[1]
            try:
                n2 = speedrange[2]
            except IndexError:
                n2 = n1
            if n1 < n2: # find id, iq, torque in fieldweakening range
                dn = r['n'][-1] - r['n'][-2]
                for nn in np.linspace(r['n'][-1]+dn, n2, nstab[1]):
                    w1 = 2*np.pi*nn*self.p
                    logger.debug("fieldweakening: n %g T %g i1max %g w1 %g u1 %g",
                                nn*60, tq, i1max, w1, u1max)
                    iq, id = self.iqd_imax_umax(i1max, w1, u1max,
                                                maxtorque=T > 0)
                    tq = self.torque_iqd(iq, id)
                    if (T > 0 and tq > 0) or (T < 0 and tq < 0):
                        r['id'].append(id)
                        r['iq'].append(iq)
                        r['n'].append(nn)
                        r['T'].append(tq)
                    else:
                        logger.warning("fieldweakening: n %g T %g i1max %g w1 %g u1 %g",
                                       nn*60, tq, i1max, w1, u1max)

            if n2 < n3:
                for nn in np.linspace(r['n'][-1]+dn/2, n3, nstab[2]):
                    w1 = 2*np.pi*nn*self.p
                    try:
                        iq, id, tq = self.mtpv(
                            w1, u1max, i1max, maxtorque=T > 0)
                        if not self._inrange((iq, id)):
                            break
                    except ValueError:
                        logger.warn("ValueError at speed %f", 60*nx)
                        break
                    r['id'].append(id)
                    r['iq'].append(iq)
                    r['n'].append(nn)
                    r['T'].append(tq)

        else:
            for t, nx in zip(T, n):
                w1 = 2*np.pi*nx*self.p
                iq, id, tq = self.iqd_torque_umax(t, w1, u1max)
                r['id'].append(id)
                r['iq'].append(iq)
                r['T'].append(tq)
                r['n'].append(nx)

        for nx, iq, id in zip(r['n'], r['iq'], r['id']):
            w1 = 2*np.pi*nx*self.p
            uq, ud = self.uqd(w1, iq, id)
            r['uq'].append(uq)
            r['ud'].append(ud)
            r['u1'].append(la.norm((ud, uq))/np.sqrt(2.0))
            r['i1'].append(la.norm((id, iq))/np.sqrt(2.0))
            r['beta'].append(np.arctan2(id, iq)/np.pi*180.)
            r['gamma'].append(np.arctan2(ud, uq)/np.pi*180.)

            r['phi'].append(r['beta'][-1] - r['gamma'][-1])
            r['cosphi'].append(np.cos(r['phi'][-1]/180*np.pi))

        pmech = np.array([2*np.pi*nx*tq for nx, tq in zip(r['n'], r['T'])])
        f1 = np.array(r['n'])*self.p
        plfe1 = self.iqd_plfe1(np.array(r['iq']), np.array(r['id']), f1)
        plfe2 = self.iqd_plfe2(np.array(r['iq']), np.array(r['id']), f1)
        plmag = self.iqd_plmag(np.array(r['iq']), np.array(r['id']), f1)
        plfe = plfe1 + plfe2 + plmag
        plcu = self.betai1_plcu(np.array(r['i1']), 2*np.pi*f1)
        pltotal = plfe + plcu
        r['pmech'] = pmech.tolist()
        r['plfe'] = plfe.tolist()
        r['plcu'] = plcu.tolist()
        r['losses'] = pltotal.tolist()
        if pmech.any():
            p1 = pmech + pltotal
            if np.abs(pmech[0]) < 1e-12:
                r['eta'] = [0]
                i = 1
            else:
                r['eta'] = []
                i = 0
            if np.all(abs(p1[i:]) > abs(pmech[i:])):
                r['eta'] += (pmech[i:]/(p1[i:])).tolist()
            else:
                r['eta'] += (p1[i:]/pmech[i:]).tolist()
        return r

    def i1beta_characteristics(self, n_list, i1_list, beta_list, u1max):
        """calculate i1-beta characteristics"""
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        for n, i1, beta in zip(n_list, i1_list, beta_list):
            w1 = 2*np.pi*n*self.p
            beta = beta/180*np.pi
            iq, id = iqd(beta, i1)
            uq, ud = self.uqd(w1, iq, id)
            u1 = la.norm((ud, uq))/np.sqrt(2)
            if u1 > u1max:
                logger.debug("u1 %s > %s", u1, u1max)
                beta = self.beta_u(w1, u1max, i1)
                logger.debug("beta %s", beta*180/np.pi)
                iq, id = iqd(beta, i1)
                logger.debug("beta %s id, %s iq %s", beta*180/np.pi, id, iq)
                uq, ud = self.uqd(w1, iq, id)
                u1 = la.norm((ud, uq))/np.sqrt(2)
                logger.debug("ud %s uq %s --> u1 %s", ud, uq, u1)

            tq = self.torque_iqd(iq, id)

            r['id'].append(id)
            r['iq'].append(iq)

            r['uq'].append(uq)
            r['ud'].append(ud)
            r['u1'].append(u1)
            r['i1'].append(la.norm((id, iq))/np.sqrt(2))
            r['T'].append(tq)
            r['beta'].append(np.arctan2(id, iq)/np.pi*180.)
            r['gamma'].append(np.arctan2(ud, uq)/np.pi*180.)

            r['n'].append(n)
            r['phi'].append(r['beta'][-1]-r['gamma'][-1])
            r['cosphi'].append(np.cos(r['phi'][-1]/180*np.pi))
            r['pmech'].append(w1/self.p*r['T'][-1])

        r['losses'] = self.iqd_losses(
            *iqd(np.array(beta_list)/180*np.pi,
                 np.array(i1_list)),
            np.array(n_list)*self.p).tolist()
        return r

    def _inrange(self, iqd):
        i1 = np.linalg.norm(iqd)/np.sqrt(2)
        iqmin, idmin = self.iqdmin(i1)
        iqmax, idmax = self.iqdmax(i1)
        return iqmin <= iqd[0] <= iqmax and idmin <= iqd[1] <= idmax


class PmRelMachineLdq(PmRelMachine):
    """Standard set of PM machine given by i1,beta parameters:
    p number of pole pairs
    m number of phases
    psim flux in Vs (RMS)
    ld d-inductance in
    lq q-inductance in H
    r1 stator resistance
    ls stator leakage inductance in H
    beta angle i1 vs up in degrees
    i1 current in A (RMS)

    optional keyword args:
    psid D-Flux in Vs (RMS)
    psiq Q-Flux in Vs (RMS)
    """

    def __init__(self,  m, p, psim=[], ld=[], lq=[],
                 r1=0, beta=[], i1=[], ls=0, **kwargs):

        super(self.__class__, self).__init__(m, p, r1, ls)
        self.psid = None
        self.betarange = (-np.pi, np.pi)
        self.i1range = (0, np.inf)
        if np.isscalar(ld):
            self.ld = lambda b, i: ld
            self.psim = lambda b, i: psim
            self.lq = lambda b, i: lq
            logger.debug("ld %s lq %s psim %s", ld, lq, psim)
            return

        if len(ld) == 1:
            try:
                self.io = iqd(min(beta)*np.pi/360, max(i1)/2)
            except:
                self.io = (1, -1)
            self.ld = lambda b, i: ld[0]
            self.psim = lambda b, i: psim[0]
            self.lq = lambda b, i: lq[0]
            logger.debug("ld %s lq %s psim %s", ld, lq, psim)
            return

        beta = np.asarray(beta)/180.0*np.pi
        if np.any(beta[beta > np.pi]):
            beta[beta > np.pi] = beta - 2*np.pi
        self.io = iqd((np.min(beta)+max(beta))/2, np.max(i1)/2)
        kx = ky = 3
        if len(i1) < 4:
            ky = len(i1)-1
        if len(beta) < 4:
            kx = len(beta)-1
        try:
            pfe = kwargs['losses']
            self._set_losspar(pfe)
            self._losses = {k: ip.RectBivariateSpline(
                beta, i1, np.array(pfe[k]),
                kx=kx, ky=ky).ev for k in (
                    'styoke_hyst', 'stteeth_hyst',
                    'styoke_eddy', 'stteeth_eddy',
                    'rotor_hyst', 'rotor_eddy',
                    'magnet')}
        except KeyError as e:
            logger.warning("loss map missing: %s", e)
            pass
        if 'psid' in kwargs:
            self.betarange = min(beta), max(beta)
            self.i1range = (0, np.max(i1))
            self.psid = lambda x, y: ip.RectBivariateSpline(
                beta, i1, np.sqrt(2)*np.asarray(kwargs['psid']),
                kx=kx, ky=ky).ev(x, y)
            self.psiq = lambda x, y: ip.RectBivariateSpline(
                beta, i1, np.sqrt(2)*np.asarray(kwargs['psiq']),
                kx=kx, ky=ky).ev(x, y)

            return

        if len(i1) < 4 or len(beta) < 4:
            if len(i1) == len(beta):
                self.ld = lambda x, y: ip.interp2d(beta, i1, ld.T)(x, y)
                self.psim = lambda x, y: ip.interp2d(beta, i1, psim.T)(x, y)
                self.lq = lambda x, y: ip.interp2d(beta, i1, lq.T)(x, y)
                logger.debug("interp2d beta %s i1 %s", beta, i1)
                return
            elif len(i1) == 1:
                self.ld = lambda x, y: ip.InterpolatedUnivariateSpline(
                    beta, ld, k=1)(x)
                self.psim = lambda x, y: ip.InterpolatedUnivariateSpline(
                    beta, psim, k=1)(x)
                self.lq = lambda x, y: ip.InterpolatedUnivariateSpline(
                    beta, lq, k=1)(x)
                logger.debug("interpolatedunivariatespline beta %s", beta)
                return
            if len(beta) == 1:
                self.ld = lambda x, y: ip.InterpolatedUnivariateSpline(
                    i1, ld, k=1)(y)
                self.psim = lambda x, y: ip.InterpolatedUnivariateSpline(
                    i1, ld, k=1)(y)
                self.lq = lambda x, y: ip.InterpolatedUnivariateSpline(
                    i1, lq, k=1)(y)
                logger.debug("interpolatedunivariatespline i1 %s", i1)
                return

            raise ValueError("unsupported array size {}x{}".format(
                len(beta), len(i1)))

        self.betarange = min(beta), max(beta)
        self.i1range = (0, np.max(i1))
        self.ld = lambda x, y: ip.RectBivariateSpline(
            beta, i1, np.asarray(ld)).ev(x, y)
        self.psim = lambda x, y: ip.RectBivariateSpline(
            beta, i1, np.asarray(psim)).ev(x, y)
        self.lq = lambda x, y: ip.RectBivariateSpline(
            beta, i1, np.asarray(lq)).ev(x, y)
        logger.debug("rectbivariatespline beta %s i1 %s", beta, i1)

    def psi(self, iq, id):
        """return psid, psiq of currents iq, id"""
        beta, i1 = betai1(np.asarray(iq), np.asarray(id))
        logger.debug('beta %f (%f, %f) i1 %f %f',
                     beta, self.betarange[0], self.betarange[1],
                     i1, self.i1range[1])
        if (self.betarange[0] <= beta <= self.betarange[1] and
                i1 <= 1.01*self.i1range[1]):
            if self.psid:
                return (self.psid(beta, i1), self.psiq(beta, i1))

            psid = self.ld(beta, i1)*id + np.sqrt(2)*self.psim(beta, i1)
            psiq = self.lq(beta, i1)*iq
            return (psid, psiq)

        return (np.nan, np.nan)

    def iqdmin(self, i1):
        """max iq, min id for given current"""
        if self.betarange[0] <= -np.pi/2 <= self.betarange[1]:
            return iqd(-np.pi/2, i1)
        if self.betarange[1] == 0:
            return iqd(self.betarange[0], i1)
        return iqd(self.betarange[1], i1)

    def iqdmax(self, i1):
        """max iq, min id for given current"""
        if self.betarange[1] == 0:
            return iqd(self.betarange[1], i1)
        return iqd(self.betarange[0], i1)

    def betai1_plfe1(self, beta, i1, f1):
        return np.sum([
            self._losses[k](beta, i1)*(f1/self.fo)**self.plexp[k] for
            k in ('styoke_eddy', 'styoke_hyst',
                  'stteeth_eddy', 'stteeth_hyst')], axis=0)

    def iqd_plfe1(self, iq, id, f1):
        return self.betai1_plfe1(*betai1(iq, id), f1)

    def betai1_plfe2(self, beta, i1, f1):
        return np.sum([
            self._losses[k](beta, i1)*(f1/self.fo)**self.plexp[k] for
            k in ('rotor_eddy', 'rotor_hyst',)], axis=0)

    def iqd_plfe2(self, iq, id, f1):
        return self.betai1_plfe2(*betai1(iq, id), f1)

    def betai1_plmag(self, beta, i1, f1):
        return self._losses['magnet'](beta, i1)*(f1/self.fo)**2

    def iqd_plmag(self, iq, id, f1):
        return self.betai1_plmag(*betai1(iq, id), f1)


class PmRelMachinePsidq(PmRelMachine):
    """Standard set of PM machine parameters:
    p number of pole pairs
    m number of phases

    psid d-flux (Vs Peak)
    psiq q-flux (Vs Peak)
    r1 stator resistance (Ohm)
    r1 stator leakage inductance (H)
    id q current (A, Peak)
    iq q current (A, Peak)
    """

    def __init__(self, m, p, psid, psiq, r1, id, iq, ls=0, **kwargs):
        super(self.__class__, self).__init__(m, p, r1, ls)

        if isinstance(psid, (float, int)):
            self._psid = lambda id, iq: np.array([[psid]])
            self._psiq = lambda id, iq: np.array([[psiq]])
            return

        psid = np.asarray(psid)
        psiq = np.asarray(psiq)
        id = np.asarray(id)
        iq = np.asarray(iq)
        self.idrange = (min(id), max(id))
        self.iqrange = (min(iq), max(iq))
        self.betarange = (-np.pi if min(iq) < 0 else -np.pi/2,
                          0 if max(iq) > 0 else -np.pi/2)
        self.i1range = (0, np.sqrt(2)*np.min(id))
        self.io = np.max(iq)/2, np.min(id)/2

        if np.any(psid.shape < (4, 4)):
            if psid.shape[0] > 1 and psid.shape[1] > 1:
                self._psid = ip.interp2d(iq, id, psid.T)
                self._psiq = ip.interp2d(iq, id, psiq.T)
                return
            if len(id) == 1 or psid.shape[1] == 1:
                self._psid = lambda x, y: ip.InterpolatedUnivariateSpline(
                    iq, psid)(x)
                self._psiq = lambda x, y: ip.InterpolatedUnivariateSpline(
                    iq, psiq)(x)
                return
            if len(iq) == 1 or psid.shape[0] == 1:
                self._psid = lambda x, y: ip.InterpolatedUnivariateSpline(
                    id, psid)(y)
                self._psiq = lambda x, y: ip.InterpolatedUnivariateSpline(
                    id, psiq)(y)
                return
            raise ValueError("unsupported array size {}x{}".format(
                len(psid.shape[0]), psid.shape[1]))

        self._psid = lambda x, y: ip.RectBivariateSpline(
            iq, id, psid).ev(x, y)
        self._psiq = lambda x, y: ip.RectBivariateSpline(
            iq, id, psiq).ev(x, y)
        try:
            pfe = kwargs['losses']
            self._set_losspar(pfe)
            self._losses = {k: ip.RectBivariateSpline(
                iq, id, np.array(pfe[k])).ev for k in (
                'styoke_hyst', 'stteeth_hyst',
                'styoke_eddy', 'stteeth_eddy',
                'rotor_hyst', 'rotor_eddy',
                'magnet')}
        except KeyError as e:
            logger.warning("loss map missing: %s", e)
            pass

    def psi(self, iq, id):
        return (self._psid(iq, id),
                self._psiq(iq, id))

    def iqdmin(self, i1):
        """max iq, min id for given current"""
        if self.idrange[0] < 0 and self.idrange[1] <= 0:
            idmin = -np.sqrt(2)*i1
        else:
            idmin = 0
        if self.idrange[0] <= idmin/np.sqrt(2):
            iqmin = -np.sqrt(2)*i1
            if self.iqrange[0] <= iqmin:
                return (iqmin, idmin)
            return self.iqrange[0], idmin

        beta = np.arccos(self.iqrange[0]/i1/np.sqrt(2))
        iqmin = np.sqrt(2)*i1*np.sin(beta)
        if self.iqrange[0] <= iqmin:
            return (iqmin, idmin)

        return self.iqrange[0], self.idrange[0]

    def iqdmax(self, i1):
        """max iq, max id for given current"""
        iqmax = np.sqrt(2)*i1
        if iqmax <= np.max(self.iqrange):
            if np.min(self.idrange) < 0 and np.max(self.idrange) <= 0:
                idmax = 0
            else:
                idmax = np.sqrt(2)*i1
            if idmax <= np.max(self.idrange):
                return (iqmax, idmax)
            return (iqmax, np.max(self.idrange))

        beta = np.arccos(self.iqrange[1]/iqmax)
        iqmax = np.sqrt(2)*i1*np.cos(beta)
        idmax = np.sqrt(2)*i1*np.sin(beta)
        if idmax <= np.max(self.idrange):
            return (iqmax, idmax)

        return iqmax, np.max(self.idrange)

    def iqd_plfe1(self, iq, id, f1):
        return np.sum([
            self._losses[k](iq, id)*(f1/self.fo)**self.plexp[k] for
            k in ('styoke_eddy', 'styoke_hyst',
                  'stteeth_eddy', 'stteeth_hyst')], axis=0)

    def betai1_plfe1(self, beta, i1, f1):
        return self.iqd_plfe1(*iqd(beta, i1), f1)

    def iqd_plfe2(self, iq, id, f1):
        return np.sum([
            self._losses[k](iq, id)*(f1/self.fo)**self.plexp[k] for
            k in ('rotor_eddy', 'rotor_hyst',)], axis=0)

    def betai1_plfe2(self, beta, i1, f1):
        return self.iqd_plfe2(*iqd(beta, i1), f1)

    def iqd_plmag(self, iq, id, f1):
        return self._losses['magnet'](iq, id)*(f1/self.fo)**2

    def betai1_plmag(self, beta, i1, f1):
        return self.iqd_plmag(*iqd(beta, i1), f1)
