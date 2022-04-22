"""
  GT-ISE Femag API

  induction machine (electrical circuit model)

  Copyright 2021: Semafor Informatik & Energie AG, Switzerland
"""
import numpy as np
import scipy.optimize as so
import logging

EPS = 1e-13
KTH = 0.003921  # temperature coefficient of resistance
TREF = 20.0  # reference temperature of resistance

eecdefaults = dict(
    zeta1=0.2,
    zeta2=1.6,
    gam=0.7,
    kh=2,
    tcu1=20,
    rotor_mass=0,
    kfric_b=1,
    pl2v=0.5,
    tcu2=20)

logger = logging.getLogger('im')


def xiskin(w, temp, zeta):
    return zeta*np.sqrt(abs(w)/(2*np.pi)/(50*(1+KTH*(temp-TREF))))


def kskinl(xi, nl):
    if abs(xi) < EPS:
        return 1.0
    xi2 = 2*xi
    nl2 = nl*nl
    return 3 / (nl2*xi2)*(np.sinh(xi2) - np.sin(xi2)) / \
        (np.cosh(xi2)-np.cos(xi2)) + \
        ((nl2-1)/(nl2*xi)*(np.sinh(xi)+np.sin(xi)) /
            (np.cosh(xi)+np.cos(xi)))


def kskinr(xi, nl):
    if abs(xi) < EPS:
        return 1.0
    xi2 = 2*xi
    nl2 = nl*nl
    return xi*((np.sinh(xi2)+np.sin(xi2))/(np.cosh(xi2)-np.cos(xi2))) + \
        ((nl2-1) / 3 * xi2*((np.sinh(xi)-np.sin(xi)) /
                            (np.cosh(xi)+np.cos(xi))))


def resistance(r0, w, temp, zeta, gam, nh):
    xi = xiskin(w, temp, zeta)
    return r0*(1.+KTH*(temp - TREF))*(gam + kskinr(xi, nh)) / (1. + gam)


class Component:
    def __init__(self, arg):
        """initialize this object either from a json file
        or from a set of parameters"""
        parameters = arg
        if type(arg) == str:
            file = open(arg)
            parameters = json.load(file)
            file.close()
        for k in parameters.keys():
            setattr(self, k, parameters[k])

    def __str__(self):
        "return string format of this object"
        return repr(self.__dict__)

    def __repr__(self):
        "representation of this object"
        return self.__str__()


class InductionMachine(Component):
    """basic induction machine model with T equivalent circuit"""

    def __init__(self, parameters):
        Component.__init__(self, parameters)
        if 'f1type' in parameters:
            self.f1ref = self.f1type
        if 'u1type' in parameters:
            self.u1ref = self.u1type
        if hasattr(self, 'f1ref'):
            self.wref = 2*np.pi*self.f1ref
            if hasattr(self, 'u1ref'):
                self.psiref = self.u1ref/self.wref
        if 'lh' in parameters:
            self.ims = 0
            self.psiref = 1
            self.iml = 1/self.lh
            self.mexp = 0
        if 'pfe' in parameters:
            self.rh = self.m*(self.psiref*self.wref)**2/self.pfe
        for k in eecdefaults.keys():
            if not hasattr(self, k):
                setattr(self, k, eecdefaults[k])

    def rfe(self, w, psi):
        """equivalent resistance for iron losses"""
        try:
            if np.isclose(w, 0):
                return 0
            else:
                return self.m*((w*w)*(psi * psi)) / (
                    (self.fec + self.fee *
                     np.power(np.abs(psi)/self.psiref, self.fexp)) *
                    (0.75 + 0.25*(np.abs(w)/self.wref)) *
                    (np.abs(w)/self.wref) *
                    np.power((np.abs(psi)/self.psiref), 2.0))
        except AttributeError:
            pass
        return self.rh

    def imag(self, psi):
        """magnetizing current"""
        return (self.iml * np.abs(psi)/self.psiref +
                self.ims*np.power(np.abs(psi)/self.psiref, self.mexp))
        # return psi/self.lh

    def lrot(self, w):
        """rotor leakage inductivity"""
        kl2 = 1.+self.pl2v*(kskinl(
            xiskin(w, self.tcu2, self.zeta2), 1)-1)
        return self.lsigma2*kl2

    def lstat(self, w):
        """stator leakage inductivity"""
        return self.lsigma1

    def rrot(self, w):
        """rotor resistance"""
        return resistance(self.r2, w, self.tcu2, self.zeta2,
                          0.0, 1)

    def rstat(self, w):
        """stator resistance"""
        return resistance(self.r1, w, self.tcu1, self.zeta1,
                          self.gam, self.kh)
        # return self.r1

    def sigma(self, w, psi):
        """leakage coefficient"""
        lh = psi / self.imag(psi)
        return (1 - (np.power(np.abs(lh), 2.0)) /
                ((lh+self.lstat(w))*(lh+self.lrot(w))))

    def u1(self, w1, psi, wm):
        """stator voltage"""
        istat = self.i1(w1, psi, wm)
        z1 = (self.rstat(w1) + w1*self.lstat(w1)*1j)
        return w1*psi*1j + istat*z1

    def i1(self, w1, psi, wm):
        """stator current"""
        imag = self.imag(psi)
        if abs(w1) > 0:
            imag += w1*psi/self.rfe(w1, psi)
        return self.i2(w1, psi, wm) + imag

    def i2(self, w1, psi, wm):
        """rotor current"""
        w2 = w1 - self.p * wm
        if abs(w2) > EPS:
            z2 = (self.rrot(w2) + w2*self.lrot(w2)*1j)
            return w2*psi*1j/z2
        return 0

    def w1torque(self, w1, u1max, psi, wm, tload):
        """calculate difference between load and motor torque"""
        # check stator voltage
        u1 = self.u1(w1, psi, wm)
        self.psi = psi
        if abs(u1) > u1max:  # must adjust flux
            self.psi = so.bisect(
                lambda psi: u1max - abs(self.u1(w1, psi, wm)),
                psi, 0.1*psi)
        return self.torque(w1, self.psi, wm)-tload

    def pullouttorque(self, w1, u1):
        """pull out torque"""
        sk = self.sk(w1, u1/w1)
        w2 = sk*w1
        r2 = self.rrot(w2)
        psi = u1/w1
        lh = psi/self.imag(psi)
        x1 = w1*(lh+self.lstat(w1))
        x2 = w1*(lh+self.lrot(w2))
        r1 = self.rstat(w1)
        sigma = self.sigma(w1, u1/w1)
        return self.m*self.p * u1**2/w1 * (1 - sigma) / \
            ((r1**2 + x1**2)*r2/(sk * x1 * x2) +
             sk * x2*(r2**2 + sigma**2 * x1**2)/(r2 * x1) + 2*r1*(1-sigma))

    def sk(self, w1, psi):
        """pullout slip"""
        r2 = self.rrot(0.)
        lh = psi/self.imag(psi)
        x1 = w1*(lh+self.lstat(w1))
        x2 = w1*(lh+self.lrot(0.))
        r1 = self.rstat(w1)
        return r2/x2*np.sqrt((r1**2+x1**2)/(
            self.sigma(w1, psi)**2*x1**2 + r1**2))

    def torque(self, w1, psi, wm):
        """electric torque (in airgap)"""
        w2 = w1-self.p*wm
        if np.isclose(w2, 0):
            return 0.
        r2 = self.rrot(w2)
        i2 = self.i2(w1, psi, wm)
        return self.m*self.p/w2*r2*(i2*i2.conjugate()).real

    def torqueu(self, w1, u1max, wm):
        """electric torque (in airgap)"""
        if np.isclose(w1, self.p*wm):
            return 0.
        # find psi
        psi = u1max/w1
        self.psi = so.fsolve(
            lambda psi: u1max - abs(self.u1(w1, psi, wm)),
            psi)[0]
        return self.torque(w1, self.psi, wm)

    def w1(self, u1max, psi, tload, wm):
        """calculate stator frequency with given torque and speed"""
        wsync = max(wm*self.p, 0.001)
        sk = self.rrot(0.)/(wsync*(self.lstat(wsync) + self.lrot(0.0)))
        if tload < 0:
            a, b = wsync, wsync*(1-sk)
        else:
            a, b = wsync, wsync*(1+sk)
        try:
            return so.bisect(self.w1torque, a, b,
                             args=(u1max, psi, wm, tload))
        except ValueError as ex:
            logger.debug("a, b: %s",
                         (a, b, tload,
                          self.torque(a, self.psi, wm),
                          self.torque(b, self.psi, wm)))
            raise ex

    def wmfweak(self, u1max, psi, torque):
        """return lower speed limit of field weakening range"""
        wm0 = u1max/psi/self.p
        try:
            return so.fsolve(
                lambda wm: u1max - np.abs(self.u1(
                    self.w1(u1max, psi, torque, wm), psi, wm)),
                wm0)[0]
        except ValueError as ex:
            logger.error("wm0 %s, umax: %s", (wm0, 0.9*wm0),
                         (u1max, np.abs(self.u1(
                             self.w1(u1max, psi, torque, wm0), psi, wm0)),
                          np.abs(self.u1(
                              self.w1(u1max, psi, torque, wm0), psi, wm0))))
            raise ex

    def characteristics(self, T, n, u1max, nsamples=50):
        """calculate torque speed characteristics.
        return dict with list values of
        id, iq, n, T, ud, uq, u1, i1,
        beta, gamma, phi, cosphi, pmech, n_type

        Keyword arguments:
        T -- (float) the maximum torque in Nm
        n -- (float) the maximum speed in 1/s
        u1max -- (float) the maximum voltage in V rms
        nsamples -- (optional) number of speed samples
        """

        wmType = self.wmfweak(u1max, self.psiref, T)
        pmmax = wmType*T
        KPO = 0.9
        wmPullout = so.fsolve(
            lambda wx: (self.pullouttorque(self.p*wx, u1max) - abs(pmmax/wx)),
            wmType)
        wmtab0 = np.linspace(wmType, 2*wmPullout[0])
        for wm, tq in zip(wmtab0, [pmmax/wx for wx in wmtab0]):
            logger.debug("u1 %g psi %g tq %g wm %g",
                         u1max, self.psiref, tq, wm)
            try:
                w1 = self.w1(u1max, self.psiref, tq, wm)
            except ValueError:
                wpo = KPO * pmmax/tq
                break

        wmMax = 1.5*wpo
        if n:
            wmMax = 2*np.pi*n
        wmtab0 = np.linspace(wpo, wmMax)
        for wm, tq in zip(wmtab0, [pmmax/wx**2 for wx in wmtab0]):
            logger.debug("u1 %g psi %g tq %g wm %g",
                         u1max, self.psiref, tq, wm)
            try:
                w1 = self.w1(u1max, self.psiref, tq, wm)
            except ValueError:
                wmMax = wm
                break

        logger.info("wmtype %f wpo %f wmmax %f", wmType, wpo, wmMax)

        if wmType < wpo < wmMax:
            wmtab = []
            dw = 0
            wmrange = sorted([0, wmType, wpo, wmMax])
            for a, b in zip(wmrange, wmrange[1:]):
                nx = int(round(nsamples*(b-a)/wmrange[-1]))
                dw = (b-a)/nx
                wmtab += np.linspace(a+dw, b, nx).tolist()
            wmtab[0] = 0
        elif wmMax < wmType:
            wmrange = sorted([0, wmMax])
            wmtab = np.linspace(0, wmMax, nsamples).tolist()
        else:
            wmrange = sorted([0, wmType, wmMax])
            nx = int(round(nsamples*wmType/wmMax))
            dw = (wmMax-wmType)/(nsamples-nx)
            wmtab = (np.linspace(0, wmType, nx).tolist() +
                     np.linspace(wmType+dw, wmMax, nsamples-nx).tolist())

        logger.info("Speed range %s", wmrange)

        def tload2(wm):
            if wm < wmType and wm < wpo:
                return T
            if wm < wpo:
                return pmmax/wm
            return wpo*pmmax/wm**2

        r = dict(u1=[], i1=[], T=[], cosphi=[], n=[],
                 plfe1=[], plcu1=[], plcu2=[])
        T = [tload2(wx) for wx in wmtab]
        tfric = self.kfric_b*self.rotor_mass*30e-3/np.pi
        w1tab = []
        for wm, tq in zip(wmtab, T):
            logger.debug("u1 %g psi %g tq %g wm %g",
                         u1max, self.psiref, tq, wm)
            try:
                w1 = self.w1(u1max, self.psiref, tq, wm)
                w1tab.append(w1)
                u1 = self.u1(w1, self.psi, wm)
                r['u1'].append(np.abs(u1))
                i1 = self.i1(w1, self.psi, wm)
                r['i1'].append(np.abs(i1))
                r['cosphi'].append(np.cos(np.angle(u1) - np.angle(i1)))
                r['plfe1'].append(self.m*np.abs(u1)**2/self.rfe(w1, self.psi))
                i2 = self.i2(w1, self.psi, wm)
                r['plcu1'].append(self.m*np.abs(i1)**2*self.rstat(w1))
                r['plcu2'].append(self.m*np.abs(i2)**2*self.rrot(w1-self.p*wm))
                r['T'].append(tq - tfric)
                r['n'].append(wm/2/np.pi)
            except ValueError:
                break
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
    import matplotlib.pyplot as plt
    import femagtools.plot

    mpars = dict(
        u1ref=2184/np.sqrt(3),
        f1ref=53,
        m=3,
        p=2,
        rh=524.2776,
        # lh=35.88343e-3,
        r1=0.043,
        r2=0.0346,
        lsigma1=1.15e-3,
        lsigma2=0.73361e-3,
        fec=9100,
        fee=0,
        fexp=7.0,
        pfric=2544.6,
        rexp=2.4,
        iml=79.8286,
        ims=27.5346,
        mexp=6.5)

    im = InductionMachine(mpars)
    torque = 6543
    u1max = 2184/np.sqrt(3)
    r = im.characteristics(torque, 0, u1max)
    # plt.plot(np.array(r['n'])*60, r['cosphi'])
    #    plt.plot(np.array(r['n'])*60, r['plcu'])
    #    plt.plot(np.array(r['n'])*60, r['losses'])
    femagtools.plot.characteristics(r)
    plt.show()
