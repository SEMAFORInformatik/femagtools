# -*- coding: utf-8 -*-
"""
    femagtools.machine
    ~~~~~~~~~~~~~~~~~~

    Calculate operating charactistics of a pm/rel machine

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import sys
import numpy as np
import numpy.linalg as la
import scipy.optimize as so
import scipy.interpolate as ip
import logging

logger = logging.getLogger(__name__)


def i1beta(iq, id):
    return (np.arctan2(id, iq),
            la.norm(np.array((iq, id))))  # /np.sqrt(2.0)


def idq(i1, beta):
    return i1*np.array((np.sin(beta),
                        np.cos(beta)))  # *np.sqrt(2.0)


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

    def iqd_torque(self, torque):
        """minimum d-q-current for torque"""
        res = so.minimize(lambda idq: la.norm(idq), (0, 0), method='SLSQP',
                          constraints=({'type': 'eq',
                                        'fun': lambda iqd:
                                        self.torque_iqd(*iqd) - torque}))
        return res.x

    def iqd_torque_umax(self, torque, w1, u1max):
        "d-q current and torque at stator frequency and max voltage"
        iq, id = self.iqd_torque(torque)
        #sys.stderr.write( "id {} iq {} torque {}\n".format(id,iq,torque ))
        # check voltage
        if la.norm(self.uqd(w1, iq, id)) <= u1max:
            return (iq, id)
        # decrease psi (flux weakening mode)
        return so.fsolve(lambda iqd: (la.norm(self.uqd(w1, *iqd)) - u1max,
                                      self.torque_iqd(*iqd)-torque), (iq, id))
        
    def beta_u(self, w1, u, i1):
        "beta at given frequency, voltage and current"
        return so.fsolve(lambda b:
                         la.norm(self.uqd(w1, *(idq(i1, b[0])[::-1])))-u,
                         -np.pi/3)[0]

    def i1_u(self, w1, u, beta):
        "i1 at given frequency, voltage and beta angle"
        return so.fsolve(lambda i1:
                         la.norm(self.uqd(w1, *(idq(i1, beta))[::-1]))-u,
                         0)[0]
    
    def iq_u(self, w1, u, id):
        "iq at given frequency, voltage and id current"
        i0 = np.sqrt(((u/w1)**2 - self.psid(0, id)**2) / self.lq(id, 0))
        return so.fsolve(lambda iq:
                         la.norm(self.uqd(w1, iq, id))-u, i0)[0]
    
    def id_u(self, w1, u, iq):
        "id at given frequency, voltage and iq current"
        i0 = np.sqrt(((u/w1)**2 - self.psim(0, iq) -
                      self.psiq(iq, 0)**2) / self.ld(0, iq))
        return so.fsolve(lambda id:
                         la.norm(self.uqd(w1, iq, id))-u, i0)[0]

    def w1_u(self, u, id, iq):
        "w1 at given voltage and id,iq current"
        w10 = np.sqrt(self.psid(iq, id)**2 + self.psiq(iq, id)**2)/u
        return so.fsolve(lambda w1:
                         la.norm(self.uqd(w1, iq, id))-u, w10)[0]

    def iq_torque(self, torque, id):
        "q current with given torque and d-current"
        return so.fsolve(lambda iq: self.torque_iqd(iq, id)-torque,
                         self.i0[1])[0]
    
    def id_torque(self, torque, iq):
        "d current with given torque and d-current"
        i0 = -0.1
        return so.fsolve(lambda id: self.torque_iqd(iq, id)-torque, i0)[0]
    
    def i1_torque(self, torque, beta):
        "i1 current with given torque and beta"
        i0 = 0  # torque/(self.m*self.p*self.psid(beta, 0))
        return so.fsolve(lambda i1: self.torque(beta, i1)-torque, i0)[0]

    def tfric(self, wm):
        return 0.0
    
    def characteristics(self, torque_list, n_list, u1):
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        for torque, n in zip(torque_list, n_list):
            w1 = 2*np.pi*n*self.p
            iq, id = self.iqd_torque_umax(torque, w1, u1)
            #print( 'p={} w1={} tq={} id={} iq={} u1={}'.format(self.p,w1,torque,id,iq,u1) )
            r['id'].append(id)
            r['iq'].append(iq)
            uq,ud=self.uqd(w1, iq, id)
            #print( 'uq={} ud={}'.format(uq, ud) )
            r['uq'].append(uq)
            r['ud'].append(ud)
            r['u1'].append(np.sqrt(ud**2+uq**2))  # /np.sqrt(2.0))
            r['i1'].append(np.sqrt(id**2+iq**2))  # /np.sqrt(2.0))
            tq = self.torque_iqd(iq, id)
            r['T'].append(tq)
            r['beta'].append(np.arctan2(id, iq)/np.pi*180.)
            r['gamma'].append(np.arctan2(ud, uq)/np.pi*180.)

            r['n'].append(n)
            r['phi'].append(r['beta'][-1] - r['gamma'][-1])
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
    """
    def __init__(self, m, p, psim, ld, lq, r1, beta=None, i1=None):
        super(self.__class__, self).__init__(m, p, r1)
        if isinstance(ld, (float, int)):
            self.ld = lambda b, i: ld
            self.lq = lambda b, i: lq
            self.psim = lambda b, i: psim
            logger.debug("ld %s lq %s psim %s", ld, lq, psim)
            return
        if len(beta) < 2:
            if isinstance(ld[0], list):  # TODO check this
                self.ld = lambda b, i: ld[-1][-1]
                self.lq = lambda b, i: lq[-1][-1]
                self.psim = lambda b, i: psim[-1][-1]
            else:
                self.ld = lambda b, i: ld[-1]
                self.lq = lambda b, i: lq[-1]
                self.psim = lambda b, i: psim[-1]
            logger.debug("ld %s lq %s psim %s", self.ld(beta[-1], 0),
                         self.lq(beta[-1], 0), self.psim(beta[-1], 0))
            return

        beta = np.asarray(beta)*np.pi/180.
        i1 = np.asarray(i1)
        ld = np.asarray(ld)
        lq = np.asarray(lq)
        psim = np.asarray(psim)

        if len(i1) < 4 or len(beta) < 4:
            if len(i1) == len(beta):
                self._ld = ip.interp2d(beta, i1, ld.T)
                self._lq = ip.interp2d(beta, i1, lq.T)
                self._psim = ip.interp2d(beta, i1, psim.T)
                self.ld = lambda x, y: self._ld(x, y)[0]
                self.lq = lambda x, y: self._lq(x, y)[0]
                self.psim = lambda x, y: self._psim(x, y)[0]
                logger.debug("interp2d beta %s i1 %s", beta, i1)
                return
            elif len(i1) == 1:
                self._ld = ip.InterpolatedUnivariateSpline(beta, ld, k=1)
                self._lq = ip.InterpolatedUnivariateSpline(beta, lq, k=1)
                self._psim = ip.InterpolatedUnivariateSpline(beta, psim, k=1)
                self.ld = lambda x, y: self._ld(x)
                self.lq = lambda x, y: self._lq(x)
                self.psim = lambda x, y: self._psim(x)
                logger.debug("interpolatedunivariatespline beta %s", beta)
                return
            if len(beta) == 1:
                self._ld = ip.InterpolatedUnivariateSpline(i1, ld, k=1)
                self._lq = ip.InterpolatedUnivariateSpline(i1, lq, k=1)
                self._psim = ip.InterpolatedUnivariateSpline(i1, psim, k=1)
                self.ld = lambda x, y: self._ld(y)
                self.lq = lambda x, y: self._lq(y)
                self.psim = lambda x, y: self._psim(y)
                logger.debug("interpolatedunivariatespline i1 %s", i1)
                return
            
            raise ValueError("unsupported array size {}x{}".format(len(beta), len(i1)))
            
        self._ld = ip.RectBivariateSpline(beta, i1, ld)
        self._lq = ip.RectBivariateSpline(beta, i1, lq)
        self._psim = ip.RectBivariateSpline(beta, i1, psim)
        self.ld = lambda x, y: self._ld(x, y)[0][0]
        self.lq = lambda x, y: self._lq(x, y)[0][0]
        self.psim = lambda x, y: self._psim(x, y)[0][0]
        logger.debug("rectbivariatespline beta %s i1 %s", beta, i1)
        
    def z(self, w1, beta, i1):
        "machine impedance in d-q reference frame"
        return np.array(((self.r1, w1*self.ld(beta, i1)),
                         (-w1*self.lq(beta, i1), self.r1)))
    
    def uqd(self, w1, iq, id):
        "stator voltage at freq and stator current in dq components"
        beta, i1 = i1beta(iq, id)
        return self.z(w1, beta, i1). dot(np.array((iq, id))) + \
            np.array((w1 * self.psim(beta, i1), 0))
        
    def ui1beta(self, w1, i1, beta):
        "stator voltage at freq and stator current in ampl and phase"
        id, iq = idq(i1, beta)
        return self.z(w1, beta, i1).dot(np.array((iq, id))) + \
            np.array((w1 * self.psim(beta, i1), 0))
        
    def torque_iqd(self, iq, id):
        "torque at q-d-current"
        beta, i1 = i1beta(iq, id)
        psid = self.ld(beta, i1)*id + self.psim(beta, i1)
        psiq = self.lq(beta, i1)*iq
        #sys.stderr.write( 'i1 {} beta {} id {} iq {} ld {} lq {}\n'.format(
        #    i1, beta, id, iq, self.ld(beta,i1), self.lq(beta,i1) ) )
        return self.m*self.p*(psid*iq - psiq*id)

    def torque(self, beta, i1):
        "torque at i1-beta-current"
        id, iq = idq(i1, beta)
        psid = self.ld(beta, i1)*id + self.psim(beta, i1)
        psiq = self.lq(beta, i1)*iq
        return self.m*self.p*(psid*iq - psiq*id)
    
    def psid(self, beta, i1):
        "psid at i1-beta-current"
        id, iq = idq(i1, beta)
        return self.ld(beta, i1)*id + self.psim(beta, i1)

    def psiq(self, beta, i1):
        "psiq at i1-beta-current"
        id, iq = idq(i1, beta)
        return self.lq(beta, i1)*iq

    def i1beta_characteristics(self, n_list, i1_list, beta_list, u1max):
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        for n, i1, beta in zip(n_list, i1_list, beta_list):
            w1 = 2*np.pi*n*self.p
            beta = beta/180*np.pi
            id, iq = idq(i1, beta)
            uq, ud = self.uqd(w1, iq, id)
            u1 = np.sqrt(ud**2+uq**2)
            if u1 > u1max:
                logger.debug("u1 %s > %s", u1, u1max)
                beta = self.beta_u(w1, u1max, i1)
                logger.debug("beta %s", beta*180/np.pi)
                id, iq = idq(i1, beta)
                logger.debug("beta %s id, %s iq %s", beta*180/np.pi, id, iq)
                uq, ud = self.uqd(w1, iq, id)
                u1 = np.sqrt(ud**2+uq**2)
                logger.debug("ud %s uq %s --> u1 %s", ud, uq, u1)
                
            tq = self.torque(beta, i1)
            #print( 'p={} w1={} tq={} id={} iq={} u1={}'.format(self.p,w1,torque,id,iq,u1) )
            r['id'].append(id)
            r['iq'].append(iq)

            #print( 'uq={} ud={}'.format(uq, ud) )
            r['uq'].append(uq)
            r['ud'].append(ud)
            r['u1'].append(u1)  # /np.sqrt(2.0))
            r['i1'].append(np.sqrt(id**2+iq**2))  # /np.sqrt(2.0))
            tq = self.torque_iqd(iq, id)
            r['T'].append(tq)
            r['beta'].append(np.arctan2(id, iq)/np.pi*180.)
            r['gamma'].append(np.arctan2(ud, uq)/np.pi*180.)

            r['n'].append(n)
            r['phi'].append(r['beta'][-1]-r['gamma'][-1])
            r['cosphi'].append(np.cos(r['phi'][-1]/180*np.pi))
            r['pmech'].append(w1/self.p*r['T'][-1])
        return r


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
            self.psid = lambda id, iq: psid / np.sqrt(2.0)
            self.psiq = lambda id, iq: psiq / np.sqrt(2.0)
            self.i0 = 0.0  # start value for iterations
            return

        self.i0 = (id[int(len(id)/2)],
                   iq[int(len(iq)/2)])  # start value for iterations
        
        psid = np.asarray(psid) / np.sqrt(2.0)
        psiq = np.asarray(psiq) / np.sqrt(2.0)
        id = np.asarray(id) / np.sqrt(2.0)
        iq = np.asarray(iq) / np.sqrt(2.0)
        
        if len(iq) < 4 or len(id) < 4:
            self._psid = ip.interp2d(iq, id, psid.T)
            self._psiq = ip.interp2d(iq, id, psiq.T)
            self.psid = lambda id, iq: self._psid(iq, id)[0]
            self.psiq = lambda id, iq: self._psiq(iq, id)[0]

        else:
            self._psid = ip.RectBivariateSpline(iq, id, psid)
            self._psiq = ip.RectBivariateSpline(iq, id, psiq)
            self.psid = lambda id, iq: self._psid(iq, id)[0][0]
            self.psiq = lambda id, iq: self._psiq(iq, id)[0][0]
        
    def z(self, w1):
        return np.array([[self.r1, 0],
                         [0, self.r1]])
    
    def uqd(self, w1, iq, id):
        "stator voltage at freq and stator current"
        return (self.z(w1).dot(np.array([iq, id])) +
                np.array([[0, w1],
                          [-w1, 0]]).dot(np.array(
                              [[self.psiq(id, iq)],
                               [self.psid(id, iq)]]))).T[0]
        
    def torque_iqd(self, iq, id):
        "torque at q-d-current"
        #sys.stderr.write( 'id {} iq {} psid {} psiq {}\n'.format(
        #    id, iq, self.psid(id,iq), self.psiq(id,iq) ) )
        return self.m*self.p*(self.psid(id, iq)*iq - self.psiq(id, iq)*id)
    
    def ld(self, id, iq):
        "ld at d-q-current"
        return (self.psid(id, iq) - self.psim(id, iq))/id

    def lq(self, id, iq):
        "lq at d-q-current"
        return self.psiq(id, iq)/iq

    def psim(self, id, iq):
        "psim at d-q-current"
        return self.psid(0, iq)

    def torque(self, beta, i1):
        "torque at i1-beta-current"
        id, iq = idq(i1, beta)
        return self.torque_iqd(iq, id)

def main(m, op):
#if all (k in m for k in ('psim','u1')) and \
#    all (k in op for k in ('T','n')):
    r10=m.get('r1',0)
    if 'beta' in m:
        pm = PmRelMachineLdq(3, m['p'], np.array(m['psim']),
                        np.array(m['ld']),
                        np.array(m['lq']),
                        r10,
                        np.array(m['beta'])*np.pi/180., np.array(m['i1']) )
    else:
        pm = PmRelMachinePsidq(3, m['p'],
                            np.array(m['psid'])/np.sqrt(2),
                            np.array(m['psiq'])/np.sqrt(2),
                            np.array(m['torque']),
                            m['r1'],
                            np.array(m['Id'])/np.sqrt(2), np.array(m['iq'])/np.sqrt(2) )
            
    return pm.characteristics( op['T'], op['n'], m['u1'])

    return {}

if __name__ == "__main__":
    import json
    import sys
    m=json.load(sys.stdin)

    r={}
    r=main(m['machine'], m['op'])
    json.dump( r, sys.stdout )
    sys.exit(0)
