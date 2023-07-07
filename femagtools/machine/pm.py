"""
  femagtools.pm
  ~~~~~~~~~~~~~

  PM/Rel synchronous machine (SPM, IPM, RM) electrical circuit model

  Copyright 2022: Semafor Informatik & Energie AG, Switzerland
"""
import logging
import warnings
import numpy as np
import numpy.linalg as la
from .utils import iqd, betai1, skin_resistance, dqparident, KTH
import scipy.optimize as so
import scipy.interpolate as ip
from functools import partial

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
        self.tmag = 20
        self.zeta1 = 0.2
        self.gam = 0.7
        self.kh = 2
        self.kpfe = 1 # iron loss factor
        self.kfric_b = 1
        self.rotor_mass = 0
        self.kth1 = KTH

        # overwritable function: skin_resistance(r0, w, tcu, kth)
        # Arguments:
        # r0: (float) dc resistance in Ohm at 20°C
        # w: (float) frequency in rad
        # tcu: (float) conductor temperature in °C
        # kth: (float) temperature coefficient (default 3.9 e-3)
        self.skin_resistance = None

        # TODO: need this for speedranges and idq_imax_umax mtpv only
        self.check_extrapolation = True
        for k in kwargs.keys():
            setattr(self, k, kwargs[k])
        try:
            self.tfric = self.kfric_b*self.rotor_mass*30e-3/np.pi
        except AttributeError:
            self.tfric = 0

        self.plexp = {'styoke_hyst': 1.0,
                      'stteeth_hyst': 1.0,
                      'styoke_eddy': 2.0,
                      'stteeth_eddy': 2.0,
                      'rotor_hyst': 1.0,
                      'rotor_eddy': 2.0}
        def nolosses(x, y):
            return 0
        self._losses = {k: nolosses for k in (
            'styoke_hyst', 'stteeth_hyst',
            'styoke_eddy', 'stteeth_eddy',
            'rotor_hyst', 'rotor_eddy',
            'magnet')}

    def pfric(self, n):
        """friction and windage losses"""
        return 2*np.pi*n*self.tfric

    def rstat(self, w):
        """stator resistance
        """
        if self.skin_resistance is not None:
            return self.skin_resistance(self.r1, w, self.tcu1, kth=self.kth1)

        return skin_resistance(self.r1, w, self.tcu1, self.zeta1,
                               self.gam, self.kh)

    def torque_iqd(self, iq, id):
        "torque at q-d-current"
        psid, psiq = self.psi(iq, id)
        tq = self.m*self.p/2*(psid*iq - psiq*id)
        return tq

    def iqd_torque(self, torque, iqd0=0, with_mtpa=True):
        """return minimum d-q-current for torque"""
        if np.abs(torque) < 1e-2:
            return (0, 0)
        if np.isscalar(iqd0):
            i0 = self.io
        else:
            i0 = iqd0
        if with_mtpa:
            res = so.minimize(
                lambda iqd: la.norm(iqd), i0, method='SLSQP',
                constraints=({'type': 'eq',
                              'fun': lambda iqd:
                              self.torque_iqd(*iqd) - torque}))
            if res.success:
                #raise ValueError(f'Torque {torque}, io {i0}: {res.message}')
                return res.x
            def func(i1):
                return torque - self.mtpa(i1)[2]
            i1 = so.fsolve(func, res.x[0])[0]
            return self.mtpa(i1)[:2]


        def tqiq(iq):
            return torque - self.torque_iqd(float(iq), 0)
        iq = so.fsolve(tqiq, (i0[0],))[0]
        return iq, 0, self.torque_iqd(iq, 0)

    def iqd_tmech(self, torque, n, iqd0=0, with_mtpa=True):
        """return minimum d-q-current for shaft torque"""
        if np.abs(torque) < 1e-2:
            return (0, 0)
        if np.isscalar(iqd0):
            i0 = self.io
        else:
            i0 = iqd0

        if with_mtpa:
            res = so.minimize(
                lambda iqd: la.norm(iqd), i0, method='SLSQP',
                constraints=({'type': 'eq',
                              'fun': lambda iqd:
                              self.tmech_iqd(*iqd, n) - torque}))
            if res.success:
                return res.x

            #logger.warning("n: %s, torque %s: %s %s",
            #                   60*n, torque, res.message, i0)
            # try a different approach:
            #raise ValueError(
            #    f'Torque {torque:.1f} speed {60*n:.1f} {res.message}')
            def func(i1):
                return torque - self.mtpa_tmech(i1, n)[2]
            i1 = so.fsolve(func, res.x[0])[0]
            return self.mtpa_tmech(i1, n)[:2]

        def tqiq(iq):
            return torque - self.tmech_iqd(float(iq), 0, n)
        iq = so.fsolve(tqiq, (i0[0],))[0]
        return iq, 0, self.tmech_iqd(iq, 0, n)

    def tloss_iqd(self, iq, id, n):
        """return loss torque of d-q current, iron loss correction factor
        and friction windage losses"""
        if n > 1e-3:
            f1 = self.p*n
            plfe = self.kpfe * (self.iqd_plfe1(iq, id, f1) + self.iqd_plfe2(iq, id, f1))
            pmag = self.iqd_plmag(iq, id, f1)
            return (plfe + pmag + self.pfric(n))/(2*np.pi*n)
        return 0

    def tmech_iqd(self, iq, id, n):
        """return shaft torque of d-q current and speed"""
        tq = self.torque_iqd(iq, id)
        return tq - self.tloss_iqd(iq, id, n)

    def uqd(self, w1, iq, id):
        """return uq, ud of frequency w1 and d-q current"""
        psid, psiq = self.psi(iq, id)
        uqd = (self.r1*iq + w1*(self.ls*id + psid),
               self.r1*id - w1*(self.ls*iq + psiq))
        #logger.debug('beta i1 %s u1 %f', betai1(iq, id), la.norm(uqd))
        return uqd

    def w1_umax(self, u, iq, id):
        """return frequency w1 at given voltage u and id, iq current

        Keyword arguments:
        u -- the maximum voltage (RMS)
        iq, id -- the d-q currents"""
        w10 = np.sqrt(2)*u/la.norm(self.psi(iq, id))
        return so.fsolve(
            lambda w1: la.norm(self.uqd(w1, iq, id))-u*np.sqrt(2),
            w10)[0]

    def w1_imax_umax(self, i1max, u1max):
        """return frequency w1 and torque at voltage u1max and current i1max

        Keyword arguments:
        u1max -- the maximum voltage (Vrms)
        i1max -- the maximum current (Arms)"""
        iq, id, T = self.mtpa(i1max)
        n0 = u1max/np.linalg.norm(self.psi(iq, id))/2/2/np.pi/self.p
        sign = -1 if i1max > 0 else 1
        res = so.minimize(
            lambda n: sign*self.mtpa_tmech(i1max, n)[2],
            n0,
            constraints={
                'type': 'eq',
                'fun': lambda n:
                np.sqrt(2)*u1max - la.norm(
                    self.uqd(2*np.pi*n*self.p,
                             *self.mtpa_tmech(i1max, n)[:2]))})
        return 2*np.pi*res.x[0]*self.p, self.mtpa_tmech(i1max, res.x[0])[2]

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

    def i1_tmech(self, torque, beta, n):
        "return i1 current with given torque and beta"
        i1, info, ier, mesg = so.fsolve(
            lambda i1: self.tmech_iqd(*iqd(beta, i1), n)-torque,
            self.i1range[1]/2,
            full_output=True)
        if ier == 1:
            return i1
        raise ValueError("no solution found for torque {}, beta {} io {}".format(
            torque, beta, self.io))

    def i1_torque(self, torque, beta):
        "return i1 current with given torque and beta"
        i1, info, ier, mesg = so.fsolve(
            lambda i1: self.torque_iqd(*iqd(beta, i1))-torque,
            self.i1range[1]/2,
            full_output=True)
        if ier == 1:
            return i1
        raise ValueError("no solution found for torque {}, beta {} io {}".format(
            torque, beta, self.io))

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
        id0 = self.iqd_torque(torque)[1]
        return so.fsolve(
            lambda id: self.torque_iqd(np.array([iq]), id)-torque, id0)[0]

    def iqd_tmech_umax(self, torque, w1, u1max, log=0, with_mtpa=True):
        """return d-q current and shaft torque at stator frequency and max voltage
        with minimal current"""
        n = w1/2/np.pi/self.p
        if with_mtpa:
            i0 = self.iqd_tmech(torque, n)
        else:
            i1 = self.i1_tmech(torque, 0, n)
            i0 = iqd(0, i1)

        if np.linalg.norm(self.uqd(w1, *i0))/np.sqrt(2) > u1max:
            beta, i1 = betai1(*i0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                i0 = iqd(so.fsolve(lambda b: u1max - np.linalg.norm(
                    self.uqd(w1, *iqd(b, i1)))/np.sqrt(2),
                                   beta)[0], i1)

            res = so.minimize(lambda iqd: np.linalg.norm(iqd), i0, method='SLSQP',
                          constraints=(
                              {'type': 'eq',
                               'fun': lambda iqd:
                               self.tmech_iqd(*iqd, n) - torque},
                              {'type': 'ineq',
                               'fun': lambda iqd:
                               np.sqrt(2)*u1max - la.norm(self.uqd(w1, *iqd))}))
            iq, id = res.x
        else:
            iq, id = i0
        if log:
            try:
                log((iq,id))
            except:
                pass  # logger is not correct
        logger.debug("iqd_tmech_umax w1=%f torque=%f %f iq=%f id=%f u1 u1 %f %f",
                     w1, torque, self.torque_iqd(iq, id), iq, id,
                     u1max, np.linalg.norm(
                         self.uqd(w1, iq, id))/np.sqrt(2))
        return iq, id, self.tmech_iqd(iq, id, n)

    def iqd_torque_umax(self, torque, w1, u1max, log=0, with_mtpa=True):
        """return d-q current and torque at stator frequency and max voltage
        with minimal current"""
        if with_mtpa:
            i0 = self.iqd_torque(torque)
        else:
            i1 = self.i1_torque(torque, 0)
            i0 = iqd(0, i1)

        if np.linalg.norm(self.uqd(w1, *i0))/np.sqrt(2) > u1max:
            beta, i1 = betai1(*i0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                i0 = iqd(so.fsolve(lambda b: u1max - np.linalg.norm(
                    self.uqd(w1, *iqd(b, i1)))/np.sqrt(2),
                                   beta)[0], i1)

        res = so.minimize(lambda iqd: la.norm(iqd), i0, method='SLSQP',
                          constraints=(
                              {'type': 'eq',
                               'fun': lambda iqd:
                               self.torque_iqd(*iqd) - torque},
                              {'type': 'ineq',
                               'fun': lambda iqd:
                               np.sqrt(2)*u1max - la.norm(self.uqd(w1, *iqd))}))

        if log:
            log(res.x)
        logger.debug("iqd_torque_umax w1=%f torque=%f %f iq=%f id=%f u1 u1 %f %f",
                     w1, torque, self.torque_iqd(*res.x), res.x[0], res.x[1],
                     u1max, np.linalg.norm(
                         self.uqd(w1, *res.x))/np.sqrt(2))
        return res.x[0], res.x[1], self.torque_iqd(*res.x)

    def iqd_pmech_imax_umax(self, n, P, i1max, u1max, with_mtpa, with_tmech, log=0):
        """return d-q current and shaft torque at speed n, P const and max voltage"""
        T = P / n / 2 / np.pi
        w1 = 2*np.pi *n * self.p
        logger.debug("field weakening mode %.2f kW @ %.0f rpm %.1f Nm; "
                     "u1=%.0f V; plfric=%.2f W",
                     P/1000, n*60, T, u1max, self.pfric(n))
        iq, id = self.iqd_torque_imax_umax(T, n, u1max, with_tmech=with_tmech)[:2]

        if with_tmech:
            tcon = {'type': 'eq',
                    'fun': lambda iqd:
                    self.tmech_iqd(*iqd, n) - T}
        else:
            tcon = {'type': 'eq',
                    'fun': lambda iqd:
                    self.torque_iqd(*iqd) - T}

        res = so.minimize(lambda iqd: np.linalg.norm(iqd), (iq, id),
                          method='SLSQP',
                          constraints=[tcon,
                                       {'type': 'ineq',
                                        'fun': lambda iqd:
                                        (u1max * np.sqrt(2)) - np.linalg.norm(
                                            self.uqd(w1, *iqd))}])
        if res.success:
            beta, i1 = betai1(*res.x)
            logger.debug("pconst %s i1 %.2f", res.x, betai1(*res.x)[1])
            if log:
                log(res.x)
            if i1 > abs(i1max):
                return self.iqd_imax_umax(i1max, w1, u1max, T,
                                          with_mtpv=False,
                                          with_tmech=with_tmech)
            if with_tmech:
                return *res.x, self.tmech_iqd(*res.x, n)
            else:
                return *res.x, self.torque_iqd(*res.x)

        if with_tmech:
            iq, id, tq = self.mtpv_tmech(w1, u1max, iqd0=res.x, maxtorque=T>0)
        else:
            iq, id, tq = self.mtpv(w1, u1max, iqd0=res.x, maxtorque=T>0)
        logger.debug("mtpa %s i1 %.2f", res.x, betai1(*res.x)[1])
        if log:
            log((iq, id, tq))
        return iq, id, tq

    def iqd_torque_imax_umax(self, torque, n, u1max, with_tmech=False, log=0):
        """return d-q current and torque at stator frequency w1,
        max voltage  and current"""
        if with_tmech:
            iq, id = self.iqd_tmech(torque, n)
        else:
            iq, id = self.iqd_torque(torque)
        w1 = 2*np.pi*n*self.p
        # Constant torque range
        if np.linalg.norm(self.uqd(w1, iq, id)) <= u1max*np.sqrt(2):
            if log:
                log((iq, id, torque))
            return (iq, id, torque)
        # Field weaking range
        imax = betai1(iq, id)[1]
        iq, id, tq = self.iqd_imax_umax(imax, w1, u1max, torque,
                                        with_mtpv=False, with_tmech=with_tmech)
        if log:
            log((iq, id, tq))
        return iq, id, tq

    def iqd_imax_umax(self, i1max, w1, u1max, torque, with_mtpv=True, with_tmech=True):
        """return d-q current and shaft torque at stator frequency and max voltage
        and max current (for motor operation if maxtorque else generator operation)"""

        if torque > 0:
            # -pi/2 --> 0
            b0, b1 = max(-np.pi/2, self.betarange[0]), 0
            if max(self.betarange) < b1:
                raise ValueError(
                    f"invalid betarange for maxtorque>0: {self.betarange}")
        else:
            # -pi/2 --> -pi
            b0, b1 = -np.pi/2, max(-np.pi, self.betarange[0])
            if min(self.betarange) > b0:
                raise ValueError(
                    f"invalid betarange for maxtorque<0: {self.betarange}")

        deps = 1e-6
        kmax = 100

        def u1norm(b):
            return np.linalg.norm(
                self.uqd(w1, *iqd(b, abs(i1max))))/np.sqrt(2)
        if u1norm(b1) < u1max:
            # must reduce current (torque)
            logger.debug('iqd_imax_umax must reduce torque')
            if with_tmech:
                iq, id, tq = self.iqd_tmech_umax(torque, w1, u1max)
            else:
                iq, id, tq = self.iqd_torque_umax(torque, w1, u1max)
            if not with_mtpv:
                return iq, id, tq
            beta, i1 = betai1(iq, id)
        else:
            for k in range(kmax):
                bx = b0 + (b1-b0)/2
                ux = u1norm(bx)
                if ux > u1max:
                    b1 = bx
                else:
                    b0 = bx
                if abs(b1-b0) < deps:
                    break
            beta, i1 = bx, i1max
            du = ux - u1max
            logger.debug("iqd_imax_umax n %f beta %f iter=%d du=%f",
                         w1/2/np.pi/self.p*60, beta/np.pi*180, k, du)
            if abs(du) < 0.1:
                iq, id = iqd(beta, abs(i1max))
            else:
                logger.debug('oops? iqd_imax_umax one more torque reduction')
                if with_tmech:
                    iq, id = self.iqd_tmech(torque, w1/2/np.pi/self.p,
                                            iqd(beta, i1))[:2]
                else:
                    iq, id = self.iqd_torque(torque, w1/2/np.pi/self.p,
                                            iqd(beta, i1))[:2]
        if with_mtpv:
            try:
                if with_tmech:
                    return self.mtpv_tmech(w1, u1max, iqd0=(iq, id))
                else:
                    return self.mtpv(w1, u1max, iqd0=(iq, id))
            except ValueError as e:
                logger.warning(e)
        iq, id = iqd(beta, abs(i1))
        if with_tmech:
            tq = self.tmech_iqd(iq, id, w1/2/np.pi/self.p)
        else:
            tq = self.torque_iqd(iq, id)
        logger.debug("iqd_imax_umax w1=%f torque=%f %f iq=%f id=%f u1 %f %f",
                     w1, torque, tq, iq, id, u1max, np.linalg.norm(
                            self.uqd(w1, iq, id))/np.sqrt(2))
        return iq, id, tq

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

    def mtpa_tmech(self, i1, n):
        """return iq, id, shaft torque at maximum torque of current i1"""
        sign = -1 if i1 > 0 else 1
        b0 = 0 if i1 > 0 else -np.pi
        bopt, fopt, iter, funcalls, warnflag = so.fmin(
            lambda x: sign*self.tmech_iqd(*iqd(x, abs(i1)), n), b0,
            full_output=True,
            disp=0)
        iq, id = iqd(bopt[0], abs(i1))
        return [iq, id, sign*fopt]

    def mtpv(self, w1, u1, iqd0=0, maxtorque=True):
        """return d-q-current, torque for voltage and frequency
        with maximum (maxtorque=True) or minimum torque """
        sign = -1 if maxtorque else 1
        if np.isscalar(iqd0):
            i0 = (-sign*self.i1range[1]/20, -self.i1range[1]/np.sqrt(2))
        else:
            i0 = iqd0
        n = w1/2/np.pi/self.p
        res = so.minimize(lambda iqd: sign*self.torque_iqd(*iqd), i0,
                          method='SLSQP',
                          constraints={
                              'type': 'eq',
                              'fun': lambda iqd:
                              np.sqrt(2)*u1 - la.norm(self.uqd(w1, *iqd))})
        if res['success']:
            return res.x[0], res.x[1], sign*res.fun
        raise ValueError(f"mtpv w1={w1} u1={u1} maxtorque={maxtorque} res: {res['message']}")

    def mtpv_tmech(self, w1, u1, iqd0=0, maxtorque=True):
        """return d-q-current, shaft torque for voltage and frequency
        with maximum (maxtorque=True) or minimum torque """
        sign = -1 if maxtorque else 1
        if np.isscalar(iqd0):
            i0 = (-sign*self.i1range[1]/20, -self.i1range[1]/np.sqrt(2))
        else:
            i0 = iqd0
        n = w1/2/np.pi/self.p
        res = so.minimize(lambda iqd: sign*self.tmech_iqd(
            *iqd, n), i0,
                          method='SLSQP',
                          constraints={'type': 'eq',
                                       'fun': lambda iqd:
                                       np.sqrt(2)*u1 - la.norm(
                                           self.uqd(w1, *iqd))})
        if res['success']:
            return res.x[0], res.x[1], sign*res.fun
        #logger.warning("w1=%.1f u1=%.1f maxtorque=%s %s: %s", w1, u1, maxtorque, res.x, res.message)
        raise ValueError(f"mtpv_tmech w1={w1} u1={u1} maxtorque={maxtorque} res: {res['message']}")

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
                       self.betai1_plcu(i1, 2*np.pi*f)], axis=0)

    def iqd_losses(self, iq, id, f):
        return np.sum([self.iqd_plfe1(iq, id, f),
                       self.iqd_plfe2(iq, id, f),
                       self.iqd_plmag(iq, id, f),
                       self.iqd_plcu(iq, id, 2*np.pi*f)], axis=0)

    def speedranges(self, i1max, u1max, speedmax,
                    with_pmconst, with_mtpa, with_mtpv, with_tmech):
        """calculate speed range intervals:
        1. const current MTPA (u < u1max)
        2. const voltage: flux reduction / MTPA and MTPV (if enabled)
        returns list of speed limit for each interval
        calculates with friction and windage losses if with_tmech=True
        """
        if with_tmech:
            w1type, T = self.w1_imax_umax(i1max, u1max)
        else:
            iq, id, T = self.mtpa(i1max)
            w1type = self.w1_umax(u1max, iq, id)
        Pmax = w1type/self.p*T
        w1max = 2*np.pi*speedmax*self.p
        wl, wu = [w1type, min(4*w1type, w1max)]
        if with_mtpv:
            kmax = 6
            self.check_extrapolation = False
        else:
            kmax = 0
            k = kmax
        logger.debug("speedranges w1type %f wu %f", wl, wu)
        for k in range(kmax):
            wx = wl + (wu-wl)/2
            try:
                if with_pmconst:
                    iq, id = self.iqd_pmech_imax_umax(
                        wx/np.pi/2/self.p, Pmax, i1max, u1max,
                        with_mtpa, with_tmech)[:2]
                else:
                    iq, id = self.iqd_imax_umax(i1max, wx, u1max,
                                                T, with_mtpv=False)[:2]
                i1 = betai1(iq, id)[1]
                try:
                    if with_tmech:
                        def voltw1(wx):
                            return np.sqrt(2)*i1 - np.linalg.norm(
                                self.mtpv_tmech(wx, u1max, iqd0=(iq, id),
                                                maxtorque=T > 0)[:2])
                    else:
                        def voltw1(wx):
                            return np.sqrt(2)*i1 - np.linalg.norm(
                                self.mtpv(wx, u1max, iqd0=(iq, id),
                                          maxtorque=T > 0)[:2])
                    w, _, ier, _ = so.fsolve(voltw1, wx, full_output=True)
                    logger.debug("fsolve ier %d T %f w %f w1 %f", ier, T, w, wx)
                    if ier in (1, 4, 5):
                        if abs(w[0]) <= wx:
                            self.check_extrapolation = True
                            return [w/2/np.pi/self.p
                                    for w in (w1type, w[0], w1max)]  # ['MTPA', 'MTPV']
                        if w[0] > wu:
                            wl += (wu-wl)/4
                        else:
                            wl = w[0]
                except ValueError as e:
                    logger.debug("%s: wx %f wl %f wu %f", e, wx, wl, wu)
                    wl = wx
                    pass
            except ValueError as e:
                logger.warning(e)
                wu = wx

            logger.debug("%d: wx %f wl %f wu %f --> %f",
                        k, wx, wl, wu, 100*(wu-wl)/wl)

        self.check_extrapolation = True
        w1max = min(w1max, self.w1_umax(
            u1max, *iqd(-np.pi/2, abs(i1max))))
        return [w/2/np.pi/self.p for w in (w1type, w1max)]  # ['MTPA']


    def operating_point(self, T, n, u1max, with_mtpa=True):
        """
        calculate single operating point.

        return dict with values for
            f1 -- (float) stator frequency in Hz
            u1 -- (float) stator phase voltage in V rms
            i1 -- (float) stator phase current in A rms
            beta -- (float) current angle in rad
            id -- (float) D-current in A peak
            iq -- (float) Q-current in A peak
            ud -- (float) D-voltage in A peak
            uq -- (float) Q-voltgae in A peak
            p1 -- (float) electrical power in W
            pmech -- (float) mechanical power in W
            plcu1 -- (float) stator copper losses in W
            plfe1 -- (float) stator iron losses in W
            plfe2 -- (float) rotor iron losses in W
            plmag -- (float) magnet losses in W
            plfric -- (float) friction losses in W, according to Tfric
            losses -- (float) total losses in W
            cosphi -- (float) power factor for this op
            eta -- (float) efficiency for this op
            T -- (float) torque for this op in Nm, copy of argument T
            Tfric -- (float) friction torque in Nm
            n -- (float) speed for this op in 1/s, copy of argument n
            tcu1/2 -- (float) temperature of statur/rotor in deg. C

        Keyword arguments:
            T -- (float) the output torque at the shaft in Nm
            n -- (float) the speed of the machine in 1/s
            u1max -- (float) the maximum phase voltage in V rms
        """

        r = {}  # result dit

        f1 = n * self.p
        w1 = f1 * 2 * np.pi
        iq, id, tq = self.iqd_tmech_umax(T, w1, u1max, with_mtpa=with_mtpa)
        uq, ud = self.uqd(w1, iq, id)
        u1 = la.norm((ud, uq)) / np.sqrt(2.0)
        i1 = la.norm((id, iq)) / np.sqrt(2.0)
        beta = np.arctan2(id, iq)
        gamma = np.arctan2(ud, uq)
        cosphi = np.cos(beta - gamma)
        pmech = tq * n * 2 * np.pi
        plfe1 = self.iqd_plfe1(iq, id, f1)
        plfe2 = self.iqd_plfe2(iq, id, f1)
        plmag = self.iqd_plmag(iq, id, f1)
        plfe = plfe1 + plfe2 + plmag
        plfric = self.pfric(n)
        plcu = self.betai1_plcu(i1, 2 * np.pi * f1)
        pltotal = plfe + plcu + plfric
        p1 = pmech + pltotal
        if np.abs(pmech) < 1e-12:
            eta = 0  # power to low for eta calculation
        elif p1 > pmech:
            eta = pmech/p1  # motor
        else:
            eta = p1/pmech  # generator
        r['i1'] = float(i1)
        r['u1'] = float(u1)
        r['iq'] = float(iq)
        r['id'] = float(id)
        r['beta'] = float(beta/ np.pi * 180)
        r['uq'] = float(uq)
        r['ud'] = float(ud)
        r['pmech'] = float(pmech)
        r['p1'] = float(p1)
        r['plfe1']= float(plfe1)
        r['plfe2'] = float(plfe2)
        r['plmag'] = float(plmag)
        r['plfric'] = float(plfric)
        r['losses'] = float(pltotal)
        r['T'] = float(tq)
        r['Tfric'] = float(self.tfric)
        r['n'] = float(n)
        r['f1'] = float(f1)
        r['eta'] = eta
        r['cosphi'] = cosphi
        r['t_cu1'] = self.tcu1
        r['t_mag'] = self.tmag

        return r

    def characteristics(self, T, n, u1max, nsamples=60,
                        with_mtpv=True, with_mtpa=True,
                        with_pmconst=True, with_tmech=True):
        """calculate torque speed characteristics.
        return dict with list values of
        id, iq, n, T, ud, uq, u1, i1,
        beta, gamma, phi, cosphi, pmech, n_type

        Keyword arguments:
        T -- the maximum torque or the list of torque values in Nm
        n -- the maximum speed or the list of speed values in 1/s
        u1max -- the maximum voltage in V rms
        nsamples -- (optional) number of speed samples
        with_mtpv -- (optional) use mtpv if True (default)
        with_pmconst -- (optional) keep pmech const if True (default)
        with_mtpa -- (optional) use mtpa if True (default) in const speed range, set id=0 if false
        with_tmech -- (optional) use friction and windage losses if True (default)
        """
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        if np.isscalar(T):
            if with_mtpa:
                iq, id = self.iqd_torque(T)
                i1max = betai1(iq, id)[1]
                if T < 0:
                    i1max = -i1max
            else:
                i1max = self.i1_torque(T, 0)
                iq, id = iqd(0, i1max)
            if with_tmech:
                w1, Tf = self.w1_imax_umax(i1max, u1max)
            else:
                iq, id = self.iqd_torque(T)
                Tf = T
                w1 = self.w1_umax(u1max, iq, id)
            n1 = w1/2/np.pi/self.p
            r['n_type'] = n1
            nmax = n
            logger.info("Type speed %.4f n: %.4f nmax %.1f T %.1f i1max %.1f",
                        60*n1, 60*n, 60*nmax, Tf, i1max)

            n1 = min(n1, nmax)
            if n1 < nmax:
                interv = 'MTPA',  # fieldweakening range is always MTPA
                if with_mtpa:
                    speedrange = [0] + self.speedranges(
                        i1max, u1max, nmax,
                        with_pmconst, with_mtpa, with_mtpv, with_tmech)
                    if len(speedrange) > 3:
                        interv = 'MTPA', 'MTPV'
                else:
                    speedrange = [0, n1, nmax]
            else:
                interv = []
                nmax = min(nmax, self.w1_umax(
                    u1max, *iqd(-np.pi/2, abs(i1max))))/2/np.pi/self.p
                speedrange = [0, n1, nmax]

            if speedrange[-1] < speedrange[-2]:
                speedrange = speedrange[:-1]
            logger.info("Speedrange T=%g %s", Tf, speedrange)
            n3 = speedrange[-1]
            nstab = [int(nsamples*(x1-x2)/n3)
                     for x1, x2 in zip(speedrange[1:],
                                       speedrange)]
            logger.info("sample intervals %s", nstab)
            for nx in np.linspace(0, n1, nstab[0]):
                if with_tmech:
                    iq, id = self.iqd_tmech(Tf, nx, (iq, id),
                                            with_mtpa or nx == n1)[:2]
                else:
                    iq, id = self.iqd_torque(Tf, (iq, id),
                                             with_mtpa or nx == n1)[:2]
                r['id'].append(id)
                r['iq'].append(iq)
                r['n'].append(nx)
                r['T'].append(Tf)

            Pmax = 2*np.pi*n1*Tf
            for ns, nu, iv in zip(nstab[1:], speedrange[2:], interv):
                # find id, iq, torque in fieldweakening range
                if ns > 0:
                    dn = (nu - r['n'][-1])/ns
                    logger.info("RANGE %s %d: %f -- %f",
                                iv, ns, r['n'][-1] + dn, nu)
                    try:
                        for nn in np.linspace(r['n'][-1]+dn, nu, ns):
                            w1 = 2*np.pi*nn*self.p
                            logger.debug("fieldweakening: n %g T %g i1max %g w1 %g u1 %g",
                                         nn*60, Tf, i1max, w1, u1max)
                            if iv == 'MTPA':
                                if with_pmconst:
                                    iq, id, tq = self.iqd_pmech_imax_umax(
                                        nn, Pmax, i1max, u1max,
                                        with_mtpa=with_mtpa,
                                        with_tmech=with_tmech)
                                else:
                                    iq, id, tq = self.iqd_imax_umax(
                                        i1max, w1, u1max,
                                        Tf, with_tmech=with_tmech,
                                        with_mtpv=(ns == 1))
                            else:
                                if with_tmech:
                                    iq, id, tq = self.mtpv_tmech(w1, u1max,
                                                                 maxtorque=T > 0)
                                else:
                                    iq, id, tq = self.mtpv(w1, u1max,
                                                           maxtorque=T > 0)
                            if (T > 0 and tq > 0) or (T < 0 and tq < 0):
                                r['id'].append(id)
                                r['iq'].append(iq)
                                r['n'].append(nn)
                                r['T'].append(tq)
                            else:
                                logger.warning("fieldweakening: n %g T %g tq %g i1max %g w1 %g u1 %g",
                                               nn*60, T, tq, i1max, w1, u1max)
                    except ValueError as e:
                        nmax = r['n'][-1]
                        logger.warning("%s: adjusted nmax %f T %f", e, nmax, r['T'][-1])
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

        if with_tmech:
            pmech = np.array([2*np.pi*nx*tq for nx, tq in zip(r['n'], r['T'])])
        else:
            pmech = np.array([2*np.pi*nx*(tq-self.tfric) for nx, tq in zip(r['n'], r['T'])])
        f1 = np.array(r['n'])*self.p
        plfe1 = self.iqd_plfe1(np.array(r['iq']), np.array(r['id']), f1)
        plfe2 = self.iqd_plfe2(np.array(r['iq']), np.array(r['id']), f1)
        plmag = self.iqd_plmag(np.array(r['iq']), np.array(r['id']), f1)
        plfe = plfe1 + plfe2 + plmag
        plcu = self.betai1_plcu(np.array(r['i1']), 2*np.pi*f1)
        plfw = self.pfric(2*np.pi*f1)
        pltotal = plfe + plcu + plfw
        r['pmech'] = pmech.tolist()
        r['plfe'] = plfe.tolist()
        r['plcu'] = plcu.tolist()
        r['plfw'] = plfw.tolist()
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

            tq = self.tmech_iqd(iq, id)

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

        super(self.__class__, self).__init__(m, p, r1, ls, **kwargs)
        self.psid = None
        self.betarange = (-np.pi, np.pi)
        self.i1range = (0, np.inf)
        if np.isscalar(ld):
            def constval(x, b, i):
                return x
            self.ld = partial(constval, ld)
            self.psim = partial(constval, psim)
            self.lq = partial(constval, lq)
            logger.debug("ld %s lq %s psim %s", ld, lq, psim)
            return

        if len(ld) == 1:
            try:
                self.io = iqd(min(beta)*np.pi/360, max(i1)/2)
            except:
                self.io = (1, -1)
            def constval(x, b, i):
                return x
            self.ld = partial(constval, ld[0])
            self.psim = partial(constval, psim[0])
            self.lq = partial(constval, lq[0])
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
            self.psid = ip.RectBivariateSpline(
                beta, i1, np.sqrt(2)*np.asarray(kwargs['psid']),
                kx=kx, ky=ky).ev
            self.psiq = ip.RectBivariateSpline(
                beta, i1, np.sqrt(2)*np.asarray(kwargs['psiq']),
                kx=kx, ky=ky).ev
            return

        if len(i1) < 4 or len(beta) < 4:
            if len(i1) == len(beta):
                self.ld = ip.interp2d(beta, i1, ld.T)
                self.psim = ip.interp2d(beta, i1, psim.T)
                self.lq = ip.interp2d(beta, i1, lq.T)
                logger.debug("interp2d beta %s i1 %s", beta, i1)
                return
            elif len(i1) == 1:
                def interp(x, b, i):
                    return ip.InterpolatedUnivariateSpline(beta, x, k=1)(b)
                self.ld = partial(interp, ld)
                self.psim = partial(interp, psim)
                self.lq = partial(interp, lq)
                logger.debug("interpolatedunivariatespline beta %s", beta)
                return
            if len(beta) == 1:
                def interp(x, b, i):
                    return ip.InterpolatedUnivariateSpline(i1, x, k=1)(i)
                self.ld = partial(interp, ld)
                self.psim = partial(interp, psim)
                self.lq = partial(interp, lq)
                logger.debug("interpolatedunivariatespline i1 %s", i1)
                return

            raise ValueError("unsupported array size {0}x{1}".format(
                len(beta), len(i1)))

        self.betarange = min(beta), max(beta)
        self.i1range = (0, np.max(i1))
        def interp(x, b, i):
            return ip.RectBivariateSpline(beta, i1, np.asarray(x)).ev(b, i)
        self.ld = partial(interp, ld)
        self.psim = partial(interp, psim)
        self.lq = partial(interp, lq)
        logger.debug("rectbivariatespline beta %s i1 %s", beta, i1)

    def psi(self, iq, id, tol=1e-4):
        """return psid, psiq of currents iq, id"""
        beta, i1 = betai1(np.asarray(iq), np.asarray(id))
        if np.isclose(beta, np.pi, atol=1e-4):
            beta = -np.pi
        #logger.debug('beta %f (%f, %f) i1 %f %f',
        #             beta, self.betarange[0], self.betarange[1],
        #             i1, self.i1range[1])
        if self.check_extrapolation:
            if (self.betarange[0]-tol > beta or
                self.betarange[1]+tol < beta or
                i1 > 1.01*self.i1range[1]):
                return (np.nan, np.nan)
        if self.psid:
            return (self.psid(beta, i1), self.psiq(beta, i1))

        psid = self.ld(beta, i1)*id + np.sqrt(2)*self.psim(beta, i1)
        psiq = self.lq(beta, i1)*iq
        return (psid, psiq)

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
        super(self.__class__, self).__init__(m, p, r1, ls, **kwargs)

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
                def interp(x, q, d):
                    return ip.InterpolatedUnivariateSpline(iq, x)(q)
                self._psid = partial(interp, psid)
                self._psiq = partial(interp, psiq)
                return
            if len(iq) == 1 or psid.shape[0] == 1:
                def interp(x, q, d):
                    return ip.InterpolatedUnivariateSpline(id, x)(d)
                self._psid = partial(interp, psid)
                self._psiq = partial(interp, psiq)
                return
            raise ValueError("unsupported array size {}".format(
                psid.shape))

        self._psid = ip.RectBivariateSpline(iq, id, psid).ev
        self._psiq = ip.RectBivariateSpline(iq, id, psiq).ev
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
