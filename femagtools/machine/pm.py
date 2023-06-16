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
        self.kfric_b = 1
        # TODO: need this for speedranges and idq_imax_umax mtpv only
        self.check_extrapolation = True
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
        if torque > 0:
            iqd0 = self.i1range[1]/2, 0
        else:
            iqd0 = -self.i1range[1]/2, 0
        res = so.minimize(
            lambda iqd: la.norm(iqd), self.io, method='SLSQP',
            constraints=({'type': 'eq',
                          'fun': lambda iqd:
                          self.torque_iqd(*iqd) - torque}))
        if not res.success:
            raise ValueError(f'Torque {torque} out of current range')
        return res.x

    def tmech_iqd(self, iq, id, n, kpfe, pfw):
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
            self.io[1],
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

    def iqd_torque_umax(self, torque, w1, u1max, log=0):
        """return d-q current and torque at stator frequency and max voltage
        with minmal current"""
        res = so.minimize(lambda iqd: la.norm(iqd), self.io, method='SLSQP',
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

    def iqd_power_umax(self, n, P, u1, kpfe, plfric, imax=0):
        """return d-q current and torque at speed n, P const and max voltage"""
        T = P / n / 2 / np.pi
        w1 =  n * self.p * 2 * np.pi
        logging.debug("field weakening mode {:.2f}kW @ {:.0f}rpm ({:.1f}Nm; "
                      "u1={:.0f}V; plfric={:.2f}W".format(
                          P/1000, n*60, T, u1, plfric))

        con1 = {'type': 'eq',
                'fun': lambda iqd: self.tmech_iqd(*iqd, n, kpfe, plfric) - T}
        con2 = {'type': 'ineq',
                'fun': lambda iqd: (u1 * np.sqrt(2))- la.norm(self.uqd(w1, *iqd))}
        con3 = {'type': 'ineq',
                'fun': lambda iqd: imax * np.sqrt(2) - la.norm(iqd)}

        if imax > 0:
            constraints = [con1, con2, con3]
        else:
            constraints = [con1, con2]
        iqd = self.iqd_torque_imax_umax(T, n, u1)
        res = so.minimize(lambda iqd: la.norm(iqd), (iqd[2], iqd[1]), method='SLSQP',
                          options={'maxiter': 50}, tol=0.1, constraints=constraints)
        if not res.success:
            logging.warning(f"current or voltage limits exceeded at n={n*60}rpm. Calculation terminated!")
        logging.debug(f"num iterations {res.nit}")
        iq, id = res.x
        tq = self.tmech_iqd(iq, id, n, kpfe, plfric)
        if not np.isclose(T, tq, 0.5/100) and n > 0:
            logging.warning("Field weakening: Torque calc {:.2f}Nm @ {:.0f}rpm; requested {:.2f}Nm!".format(tq, n*60, T))
        return iq, id, tq

    def iqd_torque_imax_umax(self, torque, n, u1max, log=0):
        """return d-q current and torque at stator frequency w1,
        max voltage  and current"""
        iq, id = self.iqd_torque(torque)
        w1 = 2*np.pi*n*self.p
        # Constant torque range
        if np.linalg.norm(self.uqd(w1, iq, id)) <= u1max*np.sqrt(2):
            if log:
                log((iq, id, torque))
            return (iq, id, torque)
        # Field weaking range
        imax = betai1(iq, id)[1]
        iq, id, tq = self.iqd_imax_umax(imax, w1, u1max, torque, with_mtpv=False)
        if log:
            log((iq, id, tq))
        return iq, id, tq

    def iqd_maxtorque_imax_umax(self, i1max, w1, u1max):
        """return d-q current and max torque at stator frequency w1,
        max voltage  and current"""
        sign = -1 if i1max > 0 else 1
        res = so.minimize(lambda iqd: sign*self.torque_iqd(*iqd),
                          self.mtpa(i1max)[:2],
                          method='SLSQP',
                          constraints=(
                              {'type': 'ineq',
                               'fun': lambda iqd:
                               np.sqrt(2)*abs(i1max) - betai1(*iqd)[1]},
                              {'type': 'eq',
                               'fun': lambda iqd:
                               np.sqrt(2)*u1max - la.norm(self.uqd(w1, *iqd))}
                          ))
        return res.x[0], res.x[1], self.torque_iqd(*res.x)

    def iqd_imax_umax(self, i1max, w1, u1max, torque, with_mtpv=True):
        """return d-q current and torque at stator frequency and max voltage
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
            iq, id, tq = self.iqd_torque_umax(torque, w1, u1max)
            if not with_mtpv:
                return iq, id, tq
            beta, i1 = betai1(iq, id)
        else:
            for k in range(kmax):
                bx = b0 + (b1-b0)/2
                ux = u1norm(bx)
                #logger.info("%d: bx %f ux %f", k, bx, ux)
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
                iq, id = self.iqd_torque(torque)
        if with_mtpv:
            # must check mtpv
            self.check_extrapolation = False
            try:
                def voltw1(wx):
                    return np.linalg.norm(
                        self.mtpv(wx, u1max,
                                  maxtorque=torque>0)[:2]) - np.sqrt(2)*i1,
                w, _, ier, _ = so.fsolve(voltw1, w1, full_output=True)
                logger.info("3: ier %d w %f w1 %f", ier, w, w1)
                if ier in (1,5) and abs(w[0]) <= w1:
                    self.check_extrapolation = True
                    return self.mtpv(w1, u1max)
            except ValueError as e:
                logger.warning("MTPV w1=%f i1max=%f, u1max %f %s",
                               w1, i1, u1max, e)
            self.check_extrapolation = True

        iq, id = iqd(beta, abs(i1))
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

    def mtpv(self, w1, u1, iqd0=0, maxtorque=True):
        """return d-q-current, torque for voltage and frequency
        with maximum (maxtorque=True) or minimum torque """
        sign = -1 if maxtorque else 1
        if np.isscalar(iqd0):
            i0 = (-sign*self.i1range[1]/20, -self.i1range[1]/np.sqrt(2))
        else:
            i0 = iqd0
        res = so.minimize(
            lambda iqd: sign*self.torque_iqd(*iqd),
            i0, method='SLSQP',
            #bounds=((0, self.i1range[1]),
            #        (-self.i1range[1], 0)),
            constraints=(
                {'type': 'eq',
                 'fun': lambda iqd:
                 np.sqrt(2)*u1 - la.norm(self.uqd(w1, *iqd))},))
        if res['success']:
            return res.x[0], res.x[1], sign*res.fun
        raise ValueError(f"mtpv w1={w1} u1={u1} maxtorque={maxtorque} res: {res['message']}")

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

    def speedranges(self, i1max, u1max, speedmax, with_mtpv, weps=1e-4):
        """calculate speed range intervals:
        1. const current MTPA (u < u1max)
        2. const voltage: flux reduction / MTPA and MTPV (if enabled)
        returns list of speed limit for each interval
        """
        iq, id, T = self.mtpa(i1max)
        logger.debug("speedrange i1 %f T %f", i1max, T)
        w1max = 2*np.pi*speedmax*self.p
        w1type = self.w1_umax(u1max, iq, id)
        wl, wu = [w1type, w1max]
        if with_mtpv:
            kmax = 6
            self.check_extrapolation = False
        else:
            kmax = 0
            k = kmax
        for k in range(kmax):
            wx = wl + (wu-wl)/2
            try:
                iq, id = self.iqd_imax_umax(i1max, wx, u1max,
                                            T, with_mtpv=False)[:2]
                i1 = betai1(iq, id)[1]
                try:
                    def voltw1(wx):
                        return np.linalg.norm(
                            self.mtpv(wx, u1max,
                                      maxtorque=T > 0)[:2]) - np.sqrt(2)*i1
                    w, _, ier, _ = so.fsolve(voltw1, wx, full_output=True)
                    logger.debug("3: ier %d i1 %f w %f w1 %f", ier, i1, w, wx)
                    if ier in (1, 4, 5):
                        self.check_extrapolation = True
                        if abs(w[0]) <= wx:
                            return [w/2/np.pi/self.p
                                    for w in (w1type, w[0], w1max)]  # ['MTPA', 'MTPV']
                        wl = w[0]
                except ValueError as e:
                    logger.debug(e)
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


    def characteristics(self, T, n, u1max, nsamples=60,
                        with_mtpv=True, with_mtpa=True):
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
        with_mtpa -- (optional) use mtpa if True (default), disables mtpv if False
        """
        r = dict(id=[], iq=[], uq=[], ud=[], u1=[], i1=[], T=[],
                 beta=[], gamma=[], phi=[], cosphi=[], pmech=[], n=[])
        if np.isscalar(T):
            if with_mtpa:
                iq, id = self.iqd_torque(T)
                i1max = betai1(iq, id)[1]
            else:
                i1max = self.i1_torque(T, 0)
                iq, id = iqd(0, i1max)
            if T < 0:
                i1max = -i1max
            w1 = self.w1_umax(u1max, iq, id)
            n1 = w1/2/np.pi/self.p
            r['n_type'] = n1
            nmax = n
            logger.info("Type speed %f n: %f nmax %f",
                        60*n1, 60*n, 60*nmax)

            n1 = min(n1, nmax)
            if n1 < nmax:
                speedrange = [0] + self.speedranges(i1max, u1max,
                                                    nmax, with_mtpv)
                if len(speedrange) > 3:
                    interv = 'MTPA', 'MTPV'
                else:
                    interv = 'MTPA',
            else:
                speedrange = [0, n1]
                nmax = n1
                interv = []

            logger.info("Speedrange T=%g %s", T, speedrange)
            n3 = speedrange[-1]
            nstab = [int(nsamples*(x1-x2)/n3)
                     for x1, x2 in zip(speedrange[1:],
                                       speedrange)]
            logger.info("sample intervals %s", nstab)
            for nx in np.linspace(0, n1, nstab[0]):
                r['id'].append(id)
                r['iq'].append(iq)
                r['n'].append(nx)
                r['T'].append(T)

            for ns, nu, iv in zip(nstab[1:], speedrange[2:], interv):
                # find id, iq, torque in fieldweakening range
                if ns == 0:
                    ns = 1
                dn = (nu - r['n'][-1])/ns
                logger.info("RANGE %s %d: %f -- %f",
                            iv, ns, r['n'][-1] + dn, nu)
                try:
                    for nn in np.linspace(r['n'][-1]+dn, nu, ns):
                        w1 = 2*np.pi*nn*self.p
                        logger.debug("fieldweakening: n %g T %g i1max %g w1 %g u1 %g",
                                     nn*60, T, i1max, w1, u1max)
                        if iv == 'MTPA':
                            iq, id, tq = self.iqd_imax_umax(i1max, w1, u1max, T,
                                                            with_mtpv=False)
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
                    logger.warning("%s: adjusted nmax %f", e, nmax)
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
        #logger.debug('beta %f (%f, %f) i1 %f %f',
        #             beta, self.betarange[0], self.betarange[1],
        #             i1, self.i1range[1])
        if self.check_extrapolation:
            if (self.betarange[0] > beta or
                self.betarange[1] < beta or
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
