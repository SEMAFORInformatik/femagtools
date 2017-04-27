import logging
import numpy as np
import numpy.linalg as la

import scipy.optimize as so
import scipy.interpolate as ip


logger = logging.getLogger(__name__)


def betai1(iq, id, meshing=True):
    if meshing:
        size = [1 if np.isscalar(x) else np.asarray(x).size
                for x in (id, iq)]
        v = (np.tile(id, size[0]),
             np.repeat(iq, size[1]))
        return (np.arctan2(v[0], v[1]),
                la.norm(v, axis=0)/np.sqrt(2.0))
    return (np.arctan2(id, iq),
            la.norm((id, iq), axis=0)/np.sqrt(2.0))
    

def iqd(beta, i1, meshing=True):
    if meshing:
        size = [1 if np.isscalar(x) else np.asarray(x).size
                for x in (beta, i1)]
        bv = np.repeat(beta, size[1])
        return np.sqrt(2.0)*np.tile(i1,
                                    size[0])*np.array([np.cos(bv),
                                                       np.sin(bv)])
    return np.sqrt(2.0)*i1*np.array([np.cos(beta),
                                     np.sin(beta)])
    

class PmRelMachine(object):
    """Abstract base class for PmRelMachines

    ::param m: number of winding phases
    ::param p: number of pole pairs
    ::param r1: stator winding resistance (in Ohm)
    """
    def __init__(self, m, p, r1):
        self.p = p
        self.m = m
        self.r1 = r1
        self.io = (0, 0)
        
    def iqd_torque(self, torque):
        """minimum d-q-current for torque"""
        res = so.minimize(lambda idq: la.norm(idq), self.io, method='SLSQP',
                          constraints=({'type': 'eq',
                                        'fun': lambda iqd:
                                        self.torque_iqd(*iqd)[0][0] - torque}))
        return res.x

    def w1_u(self, u, iq, id):
        "return frequency w1 at given voltage and id, iq current"
        w10 = np.sqrt(2)*u/la.norm(self.psi(iq, id))
        return so.fsolve(lambda w1:
                         la.norm(self.uqd(w1, iq, id))-u*np.sqrt(2), w10)[0]

    def beta_u(self, w1, u, i1):
        "beta at given frequency, voltage and current"
        return so.fsolve(lambda b:
                         la.norm(self.uqd(w1, *(iqd(i1, b[0])[::-1])))-u,
                         -np.pi/3)[0]
    
    def iqd_torque_umax(self, torque, w1, u1max):
        "d-q current and torque at stator frequency and max voltage"
        iq, id = self.iqd_torque(torque)
        # check voltage
        if la.norm(self.uqd(w1, iq, id)) <= u1max*np.sqrt(2):
            return (iq, id)
        # decrease psi (flux weakening mode)
        return so.fsolve(
            lambda iqd: (la.norm(self.uqd(w1, *iqd)) - u1max*np.sqrt(2),
                         self.torque_iqd(*iqd)[0][0] - torque),
            (iq, id))

    def characteristics(self, T, n, u1max):
        """calculate torque speed characteristics"""
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[], losses=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        for t, nx in zip(T, n):
            w1 = 2*np.pi*nx*self.p
            iq, id = self.iqd_torque_umax(t, w1, u1max)
            r['id'].append(id)
            r['iq'].append(iq)
            uq, ud = self.uqd(w1, iq, id)
            r['uq'].append(uq[0][0])
            r['ud'].append(ud[0][0])
            r['u1'].append(la.norm((ud, uq))/np.sqrt(2.0))
            r['i1'].append(la.norm((id, iq))/np.sqrt(2.0))
            tq = self.torque_iqd(iq, id)
            r['T'].append(tq)
            r['beta'].append(np.arctan2(id, iq)/np.pi*180.)
            r['gamma'].append(np.arctan2(ud[0][0], uq[0][0])/np.pi*180.)

            r['n'].append(nx)
            r['phi'].append(r['beta'][-1] - r['gamma'][-1])
            r['cosphi'].append(np.cos(r['phi'][-1]/180*np.pi))
            r['pmech'].append((2*np.pi*nx*tq)[0][0])
            #k = (cw*(p*nx/fo)**alfal + ch * (p*nx/fo)**betal)/(ch+cw)
            #r['losses'].append(k*self._losses(*betai1(iq, id))[0][0] +
            #                   la.norm((iq, id))**2/2*r1)

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
                u1 = la.norm((ud[0], uq[0]))/np.sqrt(2)
                logger.debug("ud %s uq %s --> u1 %s", ud, uq, u1)
                
            tq = self.torque_iqd(iq, id)[0][0]
            #print( 'p={} w1={} tq={} id={} iq={} u1={}'.format(self.p,w1,torque,id,iq,u1) )
            r['id'].append(id[0])
            r['iq'].append(iq[0])

            #print( 'uq={} ud={}'.format(uq, ud) )
            r['uq'].append(uq[0][0])
            r['ud'].append(ud[0][0])
            r['u1'].append(u1)
            r['i1'].append(la.norm((id, iq))/np.sqrt(2))
            r['T'].append(tq)
            r['beta'].append(np.arctan2(id[0], iq[0])/np.pi*180.)
            r['gamma'].append(np.arctan2(ud[0][0], uq[0][0])/np.pi*180.)

            r['n'].append(n)
            r['phi'].append(r['beta'][-1]-r['gamma'][-1])
            r['cosphi'].append(np.cos(r['phi'][-1]/180*np.pi))
            r['pmech'].append(w1/self.p*r['T'][-1])
        return r


class PmRelMachineLdq(PmRelMachine):
    """Standard set of PM machine given by i1,beta parameters:
    p number of pole pairs
    m number of phases
    psim flux in Vs (RMS)
    ld d-inductance in H
    lq q-inductance in H
    r1 stator resistance
    beta angle i1 vs up in degrees
    i1 current in A (RMS)

    optional keyword args:
    psid D-Flux in Vs (RMS)
    psiq Q-Flux in Vs (RMS)
    """
    def __init__(self,  m, p, psim=[], ld=[], lq=[],
                 r1=0, beta=[], i1=[], **kwargs):

        super(self.__class__, self).__init__(m, p, r1)
        if np.isscalar(ld):
            self._psid = lambda b, i: np.sqrt(2)*np.array(
                [ld*i*np.sin(b) + psim])
            self._psiq = lambda b, i: np.sqrt(2)*np.array([lq*i*np.cos(b)])
            logger.debug("ld %s lq %s psim %s", ld, lq, psim)
            return

        if len(ld) == 1:
            self.io = iqd(min(beta)*np.pi/360, max(i1)/2).ravel()
            self._psid = lambda b, i: np.sqrt(2)*np.array(
                [ld[0]*i*np.sin(b) + psim[0]])
            self._psiq = lambda b, i: np.sqrt(2)*np.array([lq[0]*i**np.cos(b)])
            logger.debug("ld %s lq %s psim %s", ld, lq, psim)
            return
        beta = np.asarray(beta)/180.0*np.pi
        self.io = iqd(np.min(beta)/2, np.max(i1)/2).ravel()
        if 'psid' not in kwargs:
            iq, id = iqd(beta, i1)
            psid = np.asarray(ld)*id.reshape((
                beta.size, len(i1))) + np.sqrt(2)*np.asarray(psim)
            psiq = np.asarray(lq)*iq.reshape((beta.size, len(i1)))
        else:
            psid = np.sqrt(2)*np.asarray(kwargs['psid'])
            psiq = np.sqrt(2)*np.asarray(kwargs['psiq'])
        if len(i1) < 4 or len(beta) < 4:
            if len(i1) == len(beta):
                self._psid = lambda x, y: np.array(
                    [ip.interp2d(beta, i1, psid.T)(x, y)])
                self._psiq = lambda x, y: np.array(
                    [ip.interp2d(beta, i1, psiq.T)(x, y)])
                logger.debug("interp2d beta %s i1 %s", beta, i1)
                return
            elif len(i1) == 1:
                self._psid = lambda x, y: np.array(
                    [ip.InterpolatedUnivariateSpline(
                        beta, psid, k=1)(x)])
                self._psiq = lambda x, y: np.array(
                    [ip.InterpolatedUnivariateSpline(
                        beta, psiq, k=1)(x)])
                logger.debug("interpolatedunivariatespline beta %s", beta)
                return
            if len(beta) == 1:
                self._psid = lambda x, y: np.array(
                    [ip.InterpolatedUnivariateSpline(
                        i1, psid, k=1)(y)])
                self._psiq = lambda x, y: np.array(
                    [ip.InterpolatedUnivariateSpline(
                        i1, psiq, k=1)(y)])
                logger.debug("interpolatedunivariatespline i1 %s", i1)
                return
            
            raise ValueError("unsupported array size {}x{}".format(
                len(beta), len(i1)))
            
        self._psid = ip.RectBivariateSpline(beta, i1, psid)
        self._psiq = ip.RectBivariateSpline(beta, i1, psiq)
        logger.debug("rectbivariatespline beta %s i1 %s", beta, i1)
    
    def torque_iqd(self, iq, id, meshing=False):
        "torque at q-d-current"
        beta, i1 = np.around(betai1(iq, id, meshing), 9)
        if meshing:
            size = [1 if np.isscalar(x) else np.asarray(x).size
                    for x in (iq, id)]
            return self.m*self.p/2*(self._psid(beta, i1)*iq.reshape(size) -
                                    self._psiq(beta, i1)*id.reshape(size))
        return self.m*self.p/2*(self._psid(beta, i1)[:, 0]*iq -
                                self._psiq(beta, i1)[:, 0]*id)

    def uqd(self, w, iq, id):
        beta, i1 = betai1(iq, id, meshing=False)
        return (self.r1*iq + w*self._psid(beta, i1),
                self.r1*id - w*self._psiq(beta, i1))

    def psi(self, iq, id):
        beta, i1 = betai1(iq, id, meshing=False)
        return (self._psid(beta, i1),
                self._psiq(beta, i1))


class PmRelMachinePsidq(PmRelMachine):
    """Standard set of PM machine parameters:
    p number of pole pairs
    m number of phases

    psid d-flux (Vs Peak)
    psiq q-flux (Vs Peak)
    r1 stator resistance
    id q current (A, Peak)
    iq q current (A, Peak)
    """

    def __init__(self, m, p, psid, psiq, r1, id, iq):
        super(self.__class__, self).__init__(m, p, r1)

        if isinstance(psid, (float, int)):
            self._psid = lambda id, iq: np.array([[psid]])
            self._psiq = lambda id, iq: np.array([[psiq]])
            return

        psid = np.asarray(psid)
        psiq = np.asarray(psiq)
        id = np.asarray(id)
        iq = np.asarray(iq)
        self.io = np.max(iq)/2, np.min(id)/2
        
        if len(iq) < 4 or len(id) < 4:
            self._psid = ip.interp2d(iq, id, psid.T)
            self._psiq = ip.interp2d(iq, id, psiq.T)

        else:
            self._psid = ip.RectBivariateSpline(iq, id, psid)
            self._psiq = ip.RectBivariateSpline(iq, id, psiq)

    def torque_iqd(self, iq, id):
        "torque at q-d-current"
        size = [1 if np.isscalar(x) else np.asarray(x).size
                for x in (iq, id)]
        return self.m*self.p/2*(self._psid(iq, id)*iq.reshape(size) -
                                self._psiq(iq, id)*id.reshape(size))

    def uqd(self, w, iq, id):
        return (self.r1*iq + w*self._psid(iq, id),
                self.r1*id - w*self._psiq(iq, id))

    def psi(self, iq, id):
        return (self._psid(iq, id),
                self._psiq(iq, id))
