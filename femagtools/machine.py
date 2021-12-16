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


def puconv(dqpar, p, NR, UR, IR):
    """convert dqpar to per unit
    arguments:
    dqpar: dict from ld-iq or psid-psiq identification
    p: pole pairs
    NR: ref speed in 1/s
    UR: ref voltage per phase in V
    IR: ref current per phase in A
    """
    WR = 2*np.pi*NR*p
    PSIR = UR/WR
    SR = 3*UR*IR
    if 'beta' in dqpar:
        dqp = dict(beta=dqpar['beta'], losses=dict())
        dqp['i1'] = np.array(dqpar['i1'])/IR
    elif 'iq' in dqpar:
        dqp = dict(iq=np.array(dqpar['iq)'])/IR*np.sqrt(2), losses=dict())
        dqp['id'] = np.array(dqpar['id'])/IR*np.sqrt(2)
    else:
        raise ValueError('invalid dqpar')
    for k in 'psid', 'psiq':
        dqp[k] = np.array(dqpar[k])/PSIR
    if 'losses' in dqpar:
        for k in ('magnet', 'styoke_hyst', 'styoke_eddy',
                  'stteeth_hyst', 'stteeth_eddy', 'rotor_hyst', 'rotor_eddy'):
            dqp['losses'][k] = np.array(dqpar['losses'][k])/SR
        dqp['losses']['speed'] = p*dqpar['losses']['speed']/WR
        dqp['losses']['ef'] = dqpar['losses']['ef']
        dqp['losses']['hf'] = dqpar['losses']['hf']
    return dqp


def __scale_losses(losses, wdg, lfe):
    if losses:
        l = {k: wdg*lfe*np.array(losses[k]) for k in (
            'styoke_hyst', 'styoke_eddy',
            'stteeth_hyst', 'stteeth_eddy',
            'rotor_hyst', 'rotor_eddy',
            'magnet')}
        l['speed'] = losses['speed']
        return l
    return {}


def dqpar_interpol(xfit, dqpars, ipkey='temperature'):
    """return interpolated parameters at temperature or exc_current

    Arguments:
      xfit -- temperature or exc_current to fit dqpars
      dqpars -- list of dict with id, iq (or i1, beta), Psid and Psiq values
      ipkey -- key (string) to interpolate
    """
    # check current range
    ckeys = (('i1', 'beta'), ('id', 'iq'))
    dqtype = 0
    fpip = {k: dqpars[0][k] for k in ckeys[dqtype]}
    fpip['losses'] = dict()
    for k in ckeys[dqtype]:
        curr = np.array([f[k] for f in dqpars], dtype=object)
        shape = curr.shape
        if curr.shape != (len(dqpars), len(curr[0])):
            raise ValueError("current range conflict")
        curr = curr.astype(float)
        if not np.array([np.allclose(curr[0], c)
                         for c in curr[1:]]).all():
            raise ValueError("current range conflict")

    try:
        speed = np.array([d['losses']['speed'] for d in dqpars])
        if (np.max(speed) - np.min(speed))/np.mean(speed) > 1e-3:
            raise ValueError("losses: speed conflict")
    except KeyError:
        pass

    sorted_dqpars = sorted(dqpars, key=lambda d: d[ipkey])
    x = [f[ipkey] for f in sorted_dqpars]
    for k in ('psid', 'psiq'):
        m = np.array([f[k] for f in sorted_dqpars]).T
        if len(x) > 2:
            fpip[k] = ip.UnivariateSpline(x, m, k=2)(xfit).T
        else:
            fpip[k] = ip.interp1d(
                x, m, fill_value='extrapolate')(xfit).T
    try:
        for k in ('styoke_hyst', 'stteeth_hyst',
                  'styoke_eddy', 'stteeth_eddy',
                  'rotor_hyst', 'rotor_eddy',
                  'magnet'):
            m = np.array([f['losses'][k] for f in sorted_dqpars]).T
            if len(x) > 2:
                fpip['losses'][k] = ip.UnivariateSpline(x, m, k=2)(xfit).T
            else:
                fpip['losses'][k] = ip.interp1d(
                    x, m, fill_value='extrapolate')(xfit).T
            fpip['losses']['speed'] = dqpars[0]['losses']['speed']
            for f in ('hf', 'ef'):
                if f in dqpars[0]['losses']:
                    fpip['losses'][f] = dqpars[0]['losses'][f]
    except KeyError:
        pass
    return x, fpip


def wdg_resistance(w1, l, d, sigma=56e3):
    """return winding resistance
    arguments:
    w1: number of turns
    l: wire length of one turn
    d: wire diameter m^2
    sigma: conductivity of wire material
    """
    a = np.pi*d**2/4
    return w1*l/sigma/a


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
            try:
                losses = __scale_losses(bch.psidq['losses'], wdg, lfe)
                losses['ef'] = bch.lossPar.get('ef', [2.0, 2.0])
                losses['eh'] = bch.lossPar.get('eh', [1.0, 1.0])
            except KeyError:
                losses = {}
            return PmRelMachinePsidq(m, p, psid, psiq, r1*lfe*wdg**2,
                                     id, iq, ls*wdg**2, losses=losses)

        if bch.type.lower().find('ld-lq-identification') >= 0:
            beta = bch.ldq['beta']
            i1 = np.array(bch.ldq['i1'])/wdg
            psid = wdg*lfe*np.array(bch.ldq['psid'])
            psiq = wdg*lfe*np.array(bch.ldq['psiq'])
            try:
                losses = __scale_losses(bch.ldq['losses'], wdg, lfe)
                losses['ef'] = bch.lossPar.get('ef', [2.0, 2.0])
                losses['eh'] = bch.lossPar.get('eh', [1.0, 1.0])
            except KeyError:
                losses = {}
            return PmRelMachineLdq(m, p, psid=psid, psiq=psiq,
                                   r1=r1*lfe*wdg**2,
                                   i1=i1, beta=beta, ls=ls*wdg**2,
                                   losses=losses)
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
        ls: leakage inductance in H
    """

    def __init__(self, m, p, r1, ls):
        self.p = p
        self.m = m
        self.r1 = r1
        self.ls = ls
        self.io = (1, -1)
        self.fo = 50.0
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
                         la.norm(self.uqd(w1, iq, id))-u*np.sqrt(2),
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
        return so.fsolve(lambda id: self.torque_iqd(iq, id)-torque, id0)[0]

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

        beta0 = max(
            self.betarange[0],
            -0.7*np.pi/2 if maxtorque else -1.4*np.pi/2)

        beta, info, ier, mesg = so.fsolve(
            lambda b: la.norm(
                self.uqd(w1, *iqd(b, i1max))) - u1max*np.sqrt(2),
            beta0,
            full_output=True)
        if ier == 1:
            return iqd(beta[0], i1max)

        return self.mtpv(w1, u1max, i1max, maxtorque)[:2]

#        raise ValueError(
#            "no solution found for imax {}, w1 {}, u1max {}".format(
#                i1max, w1, u1max))

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
        i0 = (-sign*self.i1range[1]/10, self.i1range[1]/10)

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

    def betai1_plcu(self, i1):
        return self.m*self.r1*i1**2

    def iqd_plcu(self, iq, id):
        return self.m*self.r1*(iq**2+id**2)/2

    def betai1_losses(self, beta, i1, f):
        return np.sum([self.betai1_plfe1(beta, i1, f),
                       self.betai1_plfe2(beta, i1, f),
                       self.betai1_plmag(beta, i1, f),
                       self.betai1_plcu(i1)], axis=0)

    def iqd_losses(self, iq, id, f):
        return np.sum([self.iqd_plfe1(iq, id, f),
                       self.iqd_plfe2(iq, id, f),
                       self.iqd_plmag(iq, id, f),
                       self.iqd_plcu(iq, id)], axis=0)

    def characteristics(self, T, n, u1max, nsamples=50):
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
            w1 = self.w1max(u1max, iq, id)
            w1max = self.w1max(u1max, *self.iqdmin(i1max))
            nmax = max(w1, w1max)/2/np.pi/self.p

            n1 = min(w1/2/np.pi/self.p, nmax)
            r['n_type'] = n1
            logger.info("Type speed %f n: %f nmax %f",
                        60*n1, 60*n, 60*nmax)
            try:
                w1 = self.w2_imax_umax(i1max, u1max, maxtorque=T > 0)
                n2 = w1/2/np.pi/self.p
                iqmtpv, idmtpv, tq = self.mtpv(
                    w1, u1max, i1max, maxtorque=T > 0)
                if not self._inrange((iqmtpv, idmtpv)):
                    n2 = min(nmax, n)

                logger.info("n1: %f n2: %f ",
                            60*n1, 60*n2)
            except ValueError:
                n2 = min(nmax, n)

            speedrange = sorted(
                list(set([nx for nx in [n1, n2, n] if nx <= n])))
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
            nx = n1

            for nx in np.linspace(0, n1, int(n1/dn)):
                r['id'].append(id)
                r['iq'].append(iq)
                r['n'].append(nx)
                r['T'].append(T)

            if n1 < n2:
                for nx in np.linspace(nx+dn/2, n2, int(n2/dn)):
                    w1 = 2*np.pi*nx*self.p
                    iq, id = self.iqd_imax_umax(i1max, w1, u1max,
                                                maxtorque=T > 0)
                    tq = self.torque_iqd(iq, id)
                    r['id'].append(id)
                    r['iq'].append(iq)
                    r['n'].append(nx)
                    r['T'].append(tq)
                    if T > 0 and tq < 0:
                        logger.info("2: n %g T %g i1max %g w1 %g u1 %g",
                                    nx*60, tq, i1max, w1, u1max)
            if n2 < n3:
                for nx in np.linspace(nx+dn/2, n3, int(n3/dn)):
                    w1 = 2*np.pi*nx*self.p
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
                    r['n'].append(nx)
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
        plfe = self.iqd_losses(np.array(r['iq']), np.array(r['id']),
                               np.array(r['n'])*self.p)
        plcu = self.m*self.r1*np.array(r['i1'])**2
        pltotal = plfe + plcu
        r['pmech'] = pmech.tolist()
        r['plfe'] = plfe.tolist()
        r['plcu'] = plcu.tolist()
        r['losses'] = pltotal.tolist()
        if pmech.any():
            if np.abs(pmech[0]) < 1e-12:
                r['eta'] = [np.nan] + (
                    pmech[1:]/(pmech[1:]+pltotal[1:])).tolist()
            else:
                r['eta'] = (pmech/(pmech+pltotal)).tolist()

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
        if 'psid' in kwargs:
            kx = ky = 3
            if len(i1) < 4:
                ky = len(i1)-1
            if len(beta) < 4:
                kx = len(beta)-1
            self.betarange = min(beta), max(beta)
            self.i1range = (0, np.max(i1))
            self.psid = lambda x, y: ip.RectBivariateSpline(
                beta, i1, np.sqrt(2)*np.asarray(kwargs['psid']),
                kx=kx, ky=ky).ev(x, y)
            self.psiq = lambda x, y: ip.RectBivariateSpline(
                beta, i1, np.sqrt(2)*np.asarray(kwargs['psiq']),
                kx=kx, ky=ky).ev(x, y)

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
            except KeyError:
                logger.warning("loss map missing")
                pass
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
        except KeyError:
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
