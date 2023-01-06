"""
  femagtools.sm
  ~~~~~~~~~~~~~
  wound-rotor synchronous machine (EESM) electrical circuit model

  Copyright 2022: Semafor Informatik & Energie AG, Switzerland
"""
import numpy as np
import femagtools.bch
import scipy.optimize as so
import scipy.interpolate as ip
from .utils import skin_resistance
import logging
import warnings

EPS = 1e-13

eecdefaults = dict(
    zeta1=0.2,
    zeta2=0,
    gam=0.7,
    kh=2,
    tcu1=20,
    tcu2=20,
    rotor_mass=0,
    kfric_b=1)

logger = logging.getLogger('sm')
logging.captureWarnings(True)


def _linsampl(exc, excl, a):
    """auxiliary func for linear sampling of nonlinear sequence
    arguments:
    exc: (list) nonlinear sequence of excitation current
    excl: (list) linear sequence of excitation current
    a: (array) matrix to be resampled"""
    z = []
    b, i = a.shape[1:]
    for x in range(b):
        for y in range(i):
            zl = ip.interp1d(exc, a[:, x, y], kind='linear')
            z.append(zl(excl))
    return np.array(z).T.reshape(len(excl), b, i)


def _splinterp(beta, i1, betax, i1x, a):
    """auxiliary function to increase resolution of array a
    using a cubic spline interpolation.
    arguments:
    beta: (list) n original up-i angles in rad
    i1: (list) m original current in A
    betax: (list) nx new up-i angles in rad
    i1x: (list) mx new currents
    a: (nxm array) to be interpolated"""
    f = ip.RectBivariateSpline(
        beta, i1, np.asarray(a),
        kx=3, ky=3).ev
    X, Y = np.meshgrid(betax, i1x)
    return f(X, Y).T


def _islinear(exc):
    return np.abs(np.sum(
        np.asarray(exc)**2 -
        np.linspace(exc[0], exc[-1], len(exc))**2)) < 1e-2


def gradient_respecting_bounds(bounds, fun, eps=1e-8):
    """bounds: list of tuples (lower, upper)"""
    def gradient(x):
        logger.info(x)
        fx = fun(x)
        grad = np.zeros(len(x))
        for k in range(len(x)):
            d = np.zeros(len(x))
            d[k] = eps if x[k] + eps <= bounds[k][1] else -eps
            grad[k] = (fun(x + d) - fx) / d[k]
        return grad
    return gradient


class SynchronousMachine(object):
    def __init__(self, pars, lfe=1, wdg=1):
        for k in eecdefaults.keys():
            setattr(self, k, eecdefaults[k])
        for k in pars:
            if k not in ('ldq', ):
                setattr(self, k, pars[k])

        islinear = True
        iexc = [l['ex_current'] for l in pars['ldq']]
        i1 = [i/wdg for i in pars['ldq'][-1]['i1']]
        i1x = np.linspace(i1[0], i1[-1], 12)
        beta = [b*np.pi/180 for b in pars['ldq'][-1]['beta']]
        betax = np.linspace(beta[0], beta[-1], 20)

        if _islinear(iexc):
            psid = wdg*lfe*np.array([
                _splinterp(beta, i1, betax, i1x, l['psid'])
                for l in pars['ldq']])
            psiq = wdg*lfe*np.array([
                _splinterp(beta, i1, betax, i1x, l['psiq'])
                for l in pars['ldq']])
            exc = iexc
        else:
            islinear = False
            nsamples = 10
            iexcl = np.linspace(iexc[0], iexc[-1], nsamples)
            exc = iexcl
            psid = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                [_splinterp(beta, i1, betax, i1x, l['psid'])
                 for l in pars['ldq']]))
            psiq = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                [_splinterp(beta, i1, betax, i1x, l['psiq'])
                 for l in pars['ldq']]))

        self.psidf = ip.RegularGridInterpolator(
            (exc, betax, i1x), np.sqrt(2)*psid, bounds_error=False)
        self.psiqf = ip.RegularGridInterpolator(
            (exc, betax, i1x), np.sqrt(2)*psiq, bounds_error=False)
        i1max = np.sqrt(2)*(max(i1))
        self.bounds = [(np.cos(min(beta))*i1max, i1max),
                       (-i1max, 0),
                       (iexc[0], iexc[-1])]
        self.fo = 50
        self.plexp = {'styoke_hyst': 1.0,
                      'stteeth_hyst': 1.0,
                      'styoke_eddy': 2.0,
                      'stteeth_eddy': 2.0,
                      'rotor_hyst': 1.0,
                      'rotor_eddy': 2.0}
        keys = self.plexp.keys()
        try:
            if islinear:
                pfe = {k: np.array([l['losses'][k]
                                    for l in pars['ldq']])
                       for k in keys}
            else:
                pfe = {k: _linsampl(iexc, iexcl,
                                    np.array([l['losses'][k]
                                             for l in pars['ldq']]))
                       for k in keys}
            self._losses = {k: ip.RegularGridInterpolator(
                (exc, beta, i1), lfe*np.array(pfe[k])) for k in keys}
            self._set_losspar(pars['ldq'][0]['losses']['speed'],
                              pars['ldq'][0]['losses']['ef'],
                              pars['ldq'][0]['losses']['hf'])
        except KeyError:
            logger.warning("loss map missing")
            self._losses = {k: lambda x: 0 for k in (
                'styoke_hyst', 'stteeth_hyst',
                'styoke_eddy', 'stteeth_eddy',
                'rotor_hyst', 'rotor_eddy')}

    def _set_losspar(self, speed, ef, hf):
        self.fo = speed*self.p
        self.plexp = {'styoke_hyst': hf,
                      'stteeth_hyst': hf,
                      'styoke_eddy': ef,
                      'stteeth_eddy': ef,
                      'rotor_hyst': hf,
                      'rotor_eddy': ef}

    def rstat(self, w):
        """stator resistance"""
        return skin_resistance(self.r1, w, self.tcu1, self.zeta1,
                               self.gam, self.kh)

    def rrot(self, w):
        """rotor resistance"""
        return skin_resistance(self.r2, w, self.tcu2, self.zeta2,
                               0.0, 1)

    def psi(self, iq, id, iex):
        """return psid, psiq of currents iq, id"""
        beta = np.arctan2(id, iq)
        if beta > 0:
            beta -= 2*np.pi
        i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
        try:
            return self.psidf((iex, beta, i1)), self.psiqf((iex, beta, i1))
        except ValueError as ex:
            logger.error(iex, iq, id, beta, i1)
            raise ex

    def torque_iqd(self, iq, id, iex):
        "torque at q-d-current"
        psid, psiq = self.psi(iq, id, iex)
        return self.m*self.p/2*(psid*iq - psiq*id)

    def uqd(self, w1, iq, id, iex):
        """return uq, ud of frequency w1 and d-q current"""
        psid, psiq = self.psi(iq, id, iex)
        r1 = self.rstat(w1)
        return r1*id + w1*psid, r1*iq - w1*psiq

    def plcu1(self, iqde, w1):
        r1 = self.rstat(w1)
        return 3/2*r1*(iqde[0]**2 + iqde[1]**2)

    def plcu2(self, iqde, w1=0):
        r2 = self.rrot(0)
        return r2*iqde[2]**2

    def culoss(self, iqde, w1=0):
        r2 = self.rrot(0)
        return self.plcu1(iqde, w1) + self.plcu2(iqde, 0)

    def plfe1(self, beta, i1, iex, f1):
        return np.sum([
            self._losses[k]((iex, beta, i1))*(f1/self.fo)**self.plexp[k][0]
            for k in ('styoke_eddy', 'styoke_hyst',
                      'stteeth_eddy', 'stteeth_hyst')], axis=0)

    def iqd_plcu1(self, iq, id, f1):
        return self.plcu1((iq, id), 2*np.pi*f1)

    def iqd_plcu2(self, iq, id, iex):
        return self.plcu2((iq, id, iex))

    def iqd_plfe1(self, iq, id, iex, f1):
        beta = np.arctan2(id, iq)
        beta[beta > 0] -= 2*np.pi
        i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
        return self.plfe1(beta, i1, iex, f1)

    def iqd_plfe2(self, iq, id, f1):
        return np.zeros(np.asarray(iq).shape)

    def iqd_plmag(self, iq, id, f1):
        return np.zeros(np.asarray(iq).shape)

    def iqd_torque(self, torque, disp=False, maxiter=500):
        """return currents for torque with minimal losses"""
        if torque > 0:
            startvals = self.bounds[0][1], 0, sum(self.bounds[-1])/2
        else:
            startvals = -self.bounds[0][1], 0, sum(self.bounds[-1])/2

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = so.minimize(
                self.culoss, startvals, method='SLSQP',  # trust-constr
                bounds=self.bounds,
                #            jac=gradient_respecting_bounds(self.bounds, self.culoss),
                constraints=[
                    {'type': 'eq',
                     'fun': lambda iqd: self.torque_iqd(*iqd) - torque}],
                options={'disp': disp, 'maxiter': maxiter})
        if res['success']:
            return res.x
        logger.warning("%s: torque=%f %f, io=%s",
                       res['message'], torque, self.torque_iqd(*startvals),
                       startvals)
        raise ValueError(res['message'])

    def iqd_torque_umax(self, torque, w1, u1max, io=0,
                        disp=False, maxiter=500):
        """return currents for torque with minimal losses"""
        if io == 0:
            iqde = self.iqd_torque(torque, disp, maxiter)
            if np.linalg.norm(
                    self.uqd(w1, *iqde)) <= u1max*np.sqrt(2):
                return (*iqde, torque)
            io = iqde[0], iqde[1], self.bounds[-1][0]
        res = so.minimize(
            self.culoss, io, method='SLSQP',  # trust-constr
            bounds=self.bounds,
            options={'disp': disp, 'maxiter': maxiter},
            #            jac=gradient_respecting_bounds(self.bounds, self.culoss),
            constraints=[
                {'type': 'eq',
                 'fun': lambda iqd: self.torque_iqd(*iqd) - torque},
                {'type': 'eq',
                 'fun': lambda iqd: np.linalg.norm(
                     self.uqd(w1, *iqd)) - u1max*np.sqrt(2)}])
        if res['success']:
            return (*res.x, self.torque_iqd(*res.x))
        logger.warning("%s: w1=%f torque=%f, u1max=%f, io=%s",
                       res['message'], w1, torque, u1max, io)
        raise ValueError(res['message'])

    def w1_umax(self, u, iq, id, iex):
        """return frequency w1 at given voltage u and id, iq current

        Keyword arguments:
        u -- the maximum voltage (RMS)
        iq, id -- the d-q currents"""
        w10 = np.sqrt(2)*u/np.linalg.norm(self.psi(iq, id, iex))
        return so.fsolve(
            lambda w1: np.linalg.norm(self.uqd(w1, iq, id, iex))-u*np.sqrt(2),
            w10)[0]

    def characteristics(self, T, n, u1max, nsamples=50):
        """calculate torque speed characteristics.
        return dict with list values of
        n, T, u1, i1, beta, cosphi, pmech, n_type

        Keyword arguments:
        T -- (float) the maximum torque in Nm
        n -- (float) the maximum speed in 1/s
        u1max -- (float) the maximum voltage in V rms
        nsamples -- (optional) number of speed samples
        """
        iq, id, iex = self.iqd_torque(T)
        w1type = self.w1_umax(u1max, iq, id, iex)
        wmType = w1type/self.p
        pmax = T*wmType

        def tload(wm):
            if abs(wm*T) < abs(pmax):
                return T
            return pmax/wm

        wmtab = []
        dw = 0
        wmMax = 3.5*wmType
        if n > 0:
            wmMmax = min(wmMax, 2*np.pi*n)
        if wmType > wmMax:
            wmrange = sorted([0, wmMax])
            wmtab = np.linspace(0, wmMax, nsamples).tolist()
        else:
            wmrange = sorted([0, wmType, wmMax])
            nx = int(round(nsamples*wmType/wmMax))
            dw = (wmMax-wmType)/(nsamples-nx)
            wmtab = (np.linspace(0, wmType, nx).tolist() +
                     np.linspace(wmType+dw, wmMax, nsamples-nx).tolist())

        logger.info("Speed range %s", wmrange)
        wmtab[0] = 0

        r = dict(u1=[], i1=[], id=[], iq=[], T=[], cosphi=[], n=[],
                 beta=[], plfe1=[], plcu1=[], plcu2=[])
        T = [tload(wx) for wx in wmtab]
        tfric = self.kfric_b*self.rotor_mass*30e-3/np.pi
        w1tab = []
        for wm, tq in zip(wmtab, T):
            try:
                w1 = wm*self.p
                if w1 <= w1type:
                    iq, id, iex = self.iqd_torque(tq)
                else:
                    iq, id, iex = self.iqd_torque_umax(
                        tq, w1, u1max,
                        (0.9*iq, 0.9*id,
                         min(self.bounds[-1][0], 0.9*iex)))[:-1]
                uq, ud = self.uqd(w1, iq, id, iex)
                u1 = np.linalg.norm((uq, ud))/np.sqrt(2)
                f1 = w1/2/np.pi
                r['id'].append(id)
                r['iq'].append(iq)
                r['u1'].append(u1)
                beta = np.arctan2(id, iq)
                if beta > 0:
                    beta -= 2*np.pi
                i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
                r['i1'].append(i1)
                r['beta'].append(beta/180*np.pi)
                gamma = np.arctan2(ud, uq)
                r['cosphi'].append(np.cos(gamma - beta))
                r['plfe1'].append(self.plfe1(beta, i1, iex, f1))
                r['plcu1'].append(self.m*i1**2*self.rstat(w1))
                r['plcu2'].append(iex**2*self.rrot(0))
                r['T'].append(tq-tfric)
                r['n'].append(wm/2/np.pi)
            except ValueError as ex:
                logger.debug("wm %f T %f", wm, tq)
                break

        r['plfe'] = r['plfe1']
        r['plcu'] = (np.array(r['plcu1']) + np.array(r['plcu2'])).tolist()
        r['plfric'] = [2*np.pi*n*tfric for n in r['n']]
        r['pmech'] = [2*np.pi*n*tq for n, tq in zip(r['n'], r['T'])]
        pmech = np.array(r['pmech'])
        pltotal = (np.array(r['plfe1']) + np.array(r['plfric']) +
                   np.array(r['plcu1']) + np.array(r['plcu2']))
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


if __name__ == '__main__':
    import sys
    import json
    import femagtools.plot
    import matplotlib.pyplot as plt
    with open(sys.argv[1]) as fp:
        eecpar = json.load(fp)
    m = SynchronousMachine(eecpar)
    T = 240
    u1max = 163
    nmax = 1000
    r = m.characteristics(T, 0, u1max)
    femagtools.plot.characteristics(r, 'SynchronousMachine')
    plt.show()
