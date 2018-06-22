import logging
import numpy as np
import numpy.linalg as la
from .bch import Reader

import scipy.optimize as so
import scipy.interpolate as ip


logger = logging.getLogger(__name__)


def mesh(x, y):
    """return the combined vectors x and y
    """
    size = np.asarray(x).size, np.asarray(y).size
    return (np.repeat(x, size[1]),
            np.tile(y, size[0]))


def K(d):
    """space phasor transformation matrix
    (Inverse Park Transformation) T-1 * dq
    arguments:
      d: rotation angle

    returns transformation matrix
    """
    return np.array((
        (-np.cos(d), np.sin(d)),
        (-np.cos(d-2*np.pi/3), np.sin(d-2*np.pi/3)),
        (-np.cos(d+2*np.pi/3), np.sin(d+2*np.pi/3))))


def T(d):
    """space phasor transformation matrix
    (Park Transformation) T * abc
    arguments:
      d: rotation angle

    returns transformation matrix
    """
    return np.array((
        (-np.cos(d), -np.cos(d-2*np.pi/3), -np.cos(d+2*np.pi/3)),
        (np.sin(d), np.sin(d-2*np.pi/3), np.sin(d+2*np.pi/3))))/3*2


def invpark(a, q, d):
    """ convert a dq vector to the abc reference frame
    (inverse park transformation)

    Args:
        a: rotation angle
        d: value in direct axis
        q: value in quadrature axis
    """
    if np.isscalar(a) and np.isscalar(q) and np.isscalar(d):
        return np.dot(K(a), (q, d))
    if np.isscalar(q) and np.isscalar(d):
        return np.array([K(x).dot((q, d)) for x in a]).T
    return np.array([K(x).dot((y, z)) for x, y, z in zip(a, d, q)]).T


def betai1(iq, id):
    """return beta and amplitude of dq currents"""
    return (np.arctan2(id, iq),
            la.norm((id, iq), axis=0)/np.sqrt(2.0))


def iqd(beta, i1):
    """return qd currents of beta and amplitude"""
    return np.sqrt(2.0)*i1*np.array([np.cos(beta),
                                     np.sin(beta)])


def create(bch, r1, ls, lfe=1, wdg=1):
    """create PmRelMachine from BCH

    Arguments:
      bch: BchReader or Erg object
      r1: winding resistance
      ls: winding leakage
      lfe: scale factor length
      wdg: scale factor number of windings
"""
    m = 3
    if isinstance(bch, Reader):
        p = bch.machine['p']
        if bch.type.lower().find('psid-psiq-identification') >= 0:
            id = np.array(bch.psidq['id'])/wdg
            iq = np.array(bch.psidq['iq'])/wdg
            psid = wdg*lfe*np.array(bch.psidq['psid'])
            psiq = wdg*lfe*np.array(bch.psidq['psiq'])
            return PmRelMachinePsidq(m, p, psid, psiq, r1*lfe*wdg**2,
                                     id, iq, ls*wdg**2)

        if bch.type.lower().find('ld-lq-identification') >= 0:
            beta = bch.ldq['beta']
            i1 = np.array(bch.ldq['i1'])/wdg
            psid = wdg*lfe*np.array(bch.ldq['psid'])
            psiq = wdg*lfe*np.array(bch.ldq['psiq'])
            return PmRelMachineLdq(m, p, psid=psid, psiq=psiq,
                                   r1=r1*lfe*wdg**2,
                                   i1=i1, beta=beta, ls=ls*wdg**22)
        raise ValueError("Unsupported BCH type {}".format(bch.type))
    # must be ERG type:
    p = int(round(np.sqrt(2)*bch['M_sim'][-1][-1]/(
        m*bch['Psi_d'][-1][-1] * bch['i1'][-1])))

    return PmRelMachineLdq(m, p, r1=r1*lfe*wdg**2,
                           beta=bch['beta'], i1=np.array(bch['i1'])/wdg,
                           psid=wdg*lfe*np.array(bch['Psi_d'])/np.sqrt(2),
                           psiq=wdg*lfe*np.array(bch['Psi_q'])/np.sqrt(2),
                           ls=ls*wdg**2)


class PmRelMachine(object):
    """Abstract base class for PmRelMachines

    Args:
        m: number of winding phases
        p: number of pole pairs
        r1: stator winding resistance (in Ohm)
    """
    def __init__(self, m, p, r1, ls):
        self.p = p
        self.m = m
        self.r1 = r1
        self.ls = ls
        self.io = (1, -1)
        
    def torque_iqd(self, iq, id):
        "torque at q-d-current"
        psid, psiq = self.psi(iq, id)
        tq = self.m*self.p/2*(psid*iq - psiq*id)
        return tq
    
    def iqd_torque(self, torque):
        """return minimum d-q-current for torque"""
        res = so.minimize(lambda iqd: la.norm(iqd), self.io, method='SLSQP',
                          constraints=({'type': 'eq',
                                        'fun': lambda iqd:
                                        self.torque_iqd(*iqd) - torque}))
        return res.x

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

    def w2_imax_umax(self, imax, umax):
        """return frequency at current and voltage"""
        return so.fsolve(
            lambda x: la.norm(self.mtpv(x, umax)[:2] -
                              self.iqd_imax_umax(imax, x, umax)),
            np.sqrt(2)*umax/la.norm(self.psi(*self.io)))[0]
        
    def beta_u(self, w1, u, i1):
        "beta at given frequency, voltage and current"
        return so.fsolve(lambda b:
                         la.norm(self.uqd(w1, *(iqd(b, i1))))-u*np.sqrt(2),
                         np.arctan2(self.io[1], self.io[0]))[0]
    
    def iq_u(self, w1, u, id):
        "iq at given frequency, voltage and id current"
        return so.fsolve(lambda iq:
                         la.norm(self.uqd(w1, iq, id))-u*np.sqrt(2),
                         0)[0]
    
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
        raise ValueError("no solution found for w1 {}, u1 {}, beta {}".format(
            w1, u1, beta))
    
    def id_torque(self, torque, iq):
        "return d current with given torque and d-current"
        i0 = iqd(*self.io)[1]
        return so.fsolve(lambda id: self.torque_iqd(iq, id)-torque, i0)[0]
    
    def iqd_torque_umax(self, torque, w1, u1max):
        "return d-q current and torque at stator frequency and max voltage"
        iq, id = self.iqd_torque(torque)
        # check voltage
        if la.norm(self.uqd(w1, iq, id)) <= u1max*np.sqrt(2):
            return (iq, id)
        # decrease psi (flux weakening mode), let i1 == i1max
        return so.fsolve(
            lambda iqd: (la.norm(self.uqd(w1, *iqd)) - u1max*np.sqrt(2),
                         self.torque_iqd(*iqd) - torque),
            (iq, id))

    def iqd_imax_umax(self, i1max, w1, u1max):
        """return d-q current at stator frequency and max voltage
        and max current"""

        beta0 = self.betarange[0]
        beta1 = np.sum(self.betarange)/2
        u0 = la.norm(self.uqd(w1, *iqd(beta0, i1max)))
        u1 = la.norm(self.uqd(w1, *iqd(beta1, i1max)))
        du = (u0 - u1)
        db = (beta0 - beta1)
        beta0 = beta0 + db/du*(u1max*np.sqrt(2) - u0)
        if self.betarange[0] > beta0 or self.betarange[1] < beta0:
            beta0 = self.betarange[0]

        beta, info, ier, mesg = so.fsolve(
            lambda b: la.norm(
                self.uqd(w1, *iqd(b, i1max))) - u1max*np.sqrt(2),
            beta0,
            full_output=True)
        
        if ier == 1:
            return iqd(beta[0], i1max)
        raise ValueError(
            "no solution found for imax {}, w1 {}, u1max {}".format(
                i1max, w1, u1max))
    
    def mtpa(self, i1):
        """return iq, id, torque at maximum torque of current i1"""
        bopt, fopt, iter, funcalls, warnflag = so.fmin(
            lambda x: -self.torque_iqd(*iqd(x, i1)), 0,
            full_output=True,
            disp=0)
        iq, id = iqd(bopt[0], i1)
        return [iq, id, -fopt]
   
    def mtpv(self, w1, u1):
        """return d-q-current, torque for voltage and frequency
        with maximum torque"""
        res = so.minimize(
            lambda iqd: -self.torque_iqd(*iqd),
            (0, self.io[1]), method='SLSQP',
            constraints=(
                {'type': 'ineq',
                 'fun': lambda iqd:
                 np.sqrt(2)*u1 - la.norm(self.uqd(w1, *iqd))}))
        return res.x[0], res.x[1], -res.fun  # la.norm(self.uqd(w1, *(res.x))))

    def characteristics(self, T, n, u1max, nsamples=36):
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
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[], losses=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        if np.isscalar(T):
            iq, id = self.iqd_torque(T)
            i1max = betai1(iq, id)[1]
            w1 = self.w1max(u1max, iq, id)
            nmax = max(w1,
                       self.w1max(u1max, *self.iqdmin(i1max)))/2/np.pi/self.p

            n1 = min(w1/2/np.pi/self.p, nmax)
            r['n_type'] = n1
            logger.info("Type speed %f n: %f nmax %f",
                        60*n1, 60*n, 60*nmax)
            try:
                w1 = self.w2_imax_umax(i1max, u1max)
                n2 = self.w2_imax_umax(i1max, u1max)/2/np.pi/self.p
                iqmtpv, idmtpv, tq = self.mtpv(w1, u1max)
                if self._inrange((iqmtpv, idmtpv)):
                    n2 = w1/2/np.pi/self.p
                    n3 = max(n, n2)
                else:
                    n3 = min(nmax, n)
                    n2 = n3
                    
                logger.info("n1: %f n2: %f n3 : %f ",
                            60*n1, 60*n2, 60*n3)
            except ValueError:
                n3 = min(nmax, n)
                n2 = n3

            speedrange = sorted(
                list(set([nx for nx in [n1, n2, n3] if nx <= n3])))
            n1 = speedrange[0]
            n3 = speedrange[-1]
            if n2 > n3:
                n2 = n3
            logger.info("Speed intervals %s",
                        [60*nx for nx in speedrange])
            if len(speedrange) > 2:
                nsamples = nsamples - int(speedrange[1]/(n3/nsamples))
                dn = (n3-speedrange[1])/nsamples
            else:
                dn = n3 / nsamples
            for nx in np.linspace(0, n1, int(n1/dn)):
                r['id'].append(id)
                r['iq'].append(iq)
                r['n'].append(nx)
                r['T'].append(T)

            if nx < n3:
                for nx in np.linspace(nx+dn/2, n2, int(n2)//dn):
                    w1 = 2*np.pi*nx*self.p
                    try:
                        iq, id = self.iqd_imax_umax(i1max, w1, u1max)
                    except ValueError:
                        logger.warn("ValueError at speed %f", 60*nx)
                        break
                    r['id'].append(id)
                    r['iq'].append(iq)
                    r['n'].append(nx)
                    r['T'].append(self.torque_iqd(iq, id))
                    
            if nx < n3:
                for nx in np.linspace(nx+dn/2, n3, int(n3)//dn):
                    w1 = 2*np.pi*nx*self.p
                    try:
                        iq, id, tq = self.mtpv(w1, u1max)
                        if not self._inrange((iq, id)):
                            break
                    except ValueError:
                        logger.warn("ValueError at speed %f", 60*nx)
                        break
                    r['id'].append(id)
                    r['iq'].append(iq)
                    r['n'].append(nx)
                    r['T'].append(tq)
            logger.info("Max speed %f", 60*nx)
            
        else:
            for t, nx in zip(T, n):
                w1 = 2*np.pi*nx*self.p
                iq, id = self.iqd_torque_umax(t, w1, u1max)
                r['id'].append(id)
                r['iq'].append(iq)
                tq = self.torque_iqd(iq, id)
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

        for nx, tq in zip(r['n'], r['T']):
            r['pmech'].append((2*np.pi*nx*tq))
            
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
        self.io = iqd(np.min(beta)/2, np.max(i1)/2)
        if 'psid' in kwargs:
            kx = ky = 3
            if len(i1) < 4:
                ky = len(i1)-1
            if len(beta) < 4:
                kx = len(beta)-1
            self.betarange = min(beta), max(beta)
            self.i1range = (0, np.max(i1))
            psid = np.sqrt(2)*np.asarray(kwargs['psid'])
            psiq = np.sqrt(2)*np.asarray(kwargs['psiq'])
            self.psid = lambda x, y: ip.RectBivariateSpline(
                beta, i1, psid, kx=kx, ky=ky).ev(x, y)
            self.psiq = lambda x, y: ip.RectBivariateSpline(
                beta, i1, psiq, kx=kx, ky=ky).ev(x, y)
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
    
    def psi0(self, iq, id):
        """return psid, psiq of currents iq, id"""
        beta, i1 = betai1(np.asarray(iq), np.asarray(id))
        if self.psid:
            return (self.psid(beta, i1), self.psiq(beta, i1))

        psid = self.ld(beta, i1)*id + np.sqrt(2)*self.psim(beta, i1)
        psiq = self.lq(beta, i1)*iq
        return (psid, psiq)

    def iqdmin(self, i1):
        """max iq, min id for given current"""
        return iqd(self.betarange[0], i1)
    
    def iqdmax(self, i1):
        """max iq, min id for given current"""
        return iqd(self.betarange[1], i1)


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

    def __init__(self, m, p, psid, psiq, r1, id, iq, ls=0):
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
        if np.mean(self.idrange) > 0:
            self.betarange = (0, np.pi/2)
        else:
            self.betarange = (-np.pi/2, 0)
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

    def psi(self, iq, id):
        return (self._psid(iq, id),
                self._psiq(iq, id))

    def iqdmin(self, i1):
        """max iq, min id for given current"""
        if np.min(self.idrange) < 0 and np.max(self.idrange) <= 0:
            idmin = -np.sqrt(2)*i1
        else:
            idmin = 0
        if np.min(self.idrange) <= idmin:
            iqmin = 0
            if np.min(self.iqrange) <= iqmin:
                return (iqmin, idmin)
            return np.min(self.iqrange), idmin

        i1max = np.sqrt(2)*i1
        beta = np.arccos(self.iqrange[0]/i1max)
        iqmin = np.sqrt(2)*i1*np.sin(beta)
        if np.min(self.iqrange) <= iqmin:
            return (iqmin, idmin)

        return np.min(self.iqrange), np.min(self.idrange)

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

