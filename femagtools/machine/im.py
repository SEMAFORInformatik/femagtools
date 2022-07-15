"""
  femagtools machine

  induction machine (electrical circuit model)

  Copyright 2022: Semafor Informatik & Energie AG, Switzerland
"""
import numpy as np
import scipy.optimize as so
import logging
from .utils import resistance, leakinductance, wdg_resistance
import femagtools.windings
import femagtools.parstudy
import json
import copy
import warnings

EPS = 1e-13

eecdefaults = dict(
    zeta1=0.2,
    zeta2=2.4,
    gam=0.7,
    kh=2,
    tcu1=20,
    rotor_mass=0,
    kfric_b=1,
    pl2v=0.5,
    tcu2=20)

logger = logging.getLogger('im')
logging.captureWarnings(True)


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
        """rotor leakage inductance"""
        return leakinductance(self.lsigma2, w, self.tcu2,
                              self.zeta2, 1, self.pl2v)

    def lstat(self, w):
        """stator leakage inductance"""
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
        return (1 - lh**2 /
                ((lh+self.lstat(w))*(lh+self.lrot(0))))

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

    def w1torque(self, w1, u1max, psi, wm):
        """calculate motor torque"""
        # check stator voltage
        u1 = self.u1(w1, psi, wm)
        self.psi = psi
        if abs(u1) > u1max:  # must adjust flux
            self.psi = so.fsolve(
                lambda psix: u1max - abs(self.u1(w1, psix, wm)),
                psi)[0]
        return self.torque(w1, self.psi, wm)

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.psi = so.fsolve(
                lambda psi: u1max - abs(self.u1(w1, psi, wm)),
                psi)[0]
        return self.torque(w1, self.psi, wm)

    def w1(self, u1max, psi, tload, wm):
        """calculate stator frequency with given torque and speed"""
        wsync = max(wm*self.p, 0.001)
        sk = self.rrot(0.)/(wsync*(self.lstat(wsync) + self.lrot(0.0)))
        if tload < 0:
            b = 0.999*wsync
        else:
            b = 1.001*wsync
        logger.debug("wm %s tload %s w1torque %s",
                     wm, tload, self.w1torque(b, u1max, psi, wm))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return so.fsolve(
                lambda w1: self.w1torque(w1, u1max, psi, wm) - tload, b)[0]

    def wmfweak(self, u1max, psi, torque):
        """return lower speed limit of field weakening range"""
        wm0 = u1max/psi/self.p
        return so.fsolve(
            lambda wm: u1max - np.abs(self.u1(
                self.w1(u1max, psi, torque, wm), psi, wm)),
            wm0)[0]

    def characteristics(self, T, n, u1max, nsamples=50, kpo=0.9):
        """calculate torque speed characteristics.
        return dict with list values of
        id, iq, n, T, ud, uq, u1, i1,
        beta, gamma, phi, cosphi, pmech, n_type

        Keyword arguments:
        T -- (float) the maximum torque in Nm
        n -- (float) the maximum speed in 1/s
        u1max -- (float) the maximum phase voltage in V rms
        nsamples -- (optional) number of speed samples
        """

        wmType = self.wmfweak(u1max, self.psiref, T)
        pmmax = wmType*T
        wmPullout = so.fsolve(
            lambda wx: (kpo*self.pullouttorque(self.p *
                        wx, u1max) - abs(pmmax/wx)),
            wmType)[0]
        wmtab0 = np.linspace(wmType, 3*wmPullout)
        for wm, tq in zip(wmtab0, [pmmax/wx for wx in wmtab0]):
            logger.debug("u1 %g psi %g tq %g wm %g",
                         u1max, self.psiref, tq, wm)
            try:
                w1 = self.w1(u1max, self.psiref, tq, wm)
            except ValueError:
                wmPullout = wm
                break
        wmMax = max(1.5*wmPullout, 3*abs(pmmax/T))
        if n:
            wmMax = 2*np.pi*n
        if wmMax > wmPullout:
            wmtab0 = np.linspace(wmPullout, wmMax)
            for wm, tq in zip(wmtab0, [pmmax/wx**2
                                       for wx in wmtab0]):
                logger.debug("u1 %g psi %g tq %g wm %g",
                             u1max, self.psiref, tq, wm)
                try:
                    w1 = self.w1(u1max, self.psiref, tq, wm)
                except ValueError:
                    wmMax = wm
                    break

        logger.info("wmtype %f wpo %f wmmax %f", wmType, wmPullout, wmMax)

        if wmType < wmPullout < wmMax:
            wmrange = sorted([0, wmType, wmPullout, wmMax])
        elif wmMax < wmType:
            wmrange = sorted([0, wmMax])
        else:
            wmrange = sorted([0, wmType, wmMax])

        logger.info("Speed range %s", wmrange)
        wmlin = []
        dw = 0
        for i, nx in enumerate([round(nsamples*(w1-w0)/wmMax)
                                for w1, w0 in zip(wmrange[1:],
                                                  wmrange)]):
            if nx == 1:
                nx = 2
            if nx > 0:
                lw = np.linspace(wmrange[i]+dw, wmrange[i+1], nx)
                dw = lw[-1] - lw[-2]
                wmlin.append(lw)
        if len(wmlin) > 1:
            wmtab = np.concatenate(wmlin)
        else:
            wmtab = wmlin[0]

        def tload2(wm):
            if wm < wmType and wm < wmPullout:
                return T
            if wm < wmPullout:
                if pmmax < 0:
                    return max(T, pmmax/wm)
                return min(T, pmmax/wm)
            if pmmax < 0:
                return max(wmPullout*pmmax/wm**2, T)
            return min(wmPullout*pmmax/wm**2, T)

        r = dict(u1=[], i1=[], T=[], cosphi=[], n=[],
                 plfe1=[], plcu1=[], plcu2=[])
        T = [tload2(wx) for wx in wmtab]
        tfric = self.kfric_b*self.rotor_mass*30e-3/np.pi
        w1tab = []
        for wm, tq in zip(wmtab, T):
            logger.debug("u1 %g psi %g tq %g wm %g",
                         u1max, self.psiref, tq, wm)
#            try:
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
#            except ValueError as ex:
#                break
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


def parident(workdir, engine, f1, u1, wdgcon,
             machine, magnetizingCurves, condMat=[], **kwargs):
    """return dict of parameters of equivalent circuit for induction machines

    arguments:
    workdir --directory for intermediate files
    engine -- calculation driver (multiproc, docker, condor)
    f1 -- stator voltage frequency [Hz]
    u1 -- stator voltage (line-to-line) [V]
    wdgcon -- winding connection  ('open', 'wye', 'star', 'delta')
    machine -- dict() with machine parameters
    magnetizingCurves -- list of dict() with BH curves (or directory with MC/MCV files)
    condMat -- list of dict() with conductor material properties

    optional:
    num_u1_steps: (int) number of u1 steps for the no-load and load calculation (default 6)
    u1_logspace: (bool) uses log distributed voltage samples if true (linear otherwise) (default true)
    """
    import scipy.interpolate as ip
    import scipy.optimize as so
    CON = {'open': 0, 'wye': 1, 'star': 1, 'delta': 2}
    p = machine['poles']//2
    slip = 1e-2
    u1ph = u1
    if CON[wdgcon] == 1:
        u1ph /= np.sqrt(3)
    u1max = 1.25*u1ph
    u1min = u1ph/4
    num_u1_steps = kwargs.get('num_u1_steps', 6)
    f1tab = [f1]*num_u1_steps
    u1_logspace = kwargs.get('u1_logspace', True)
    if u1_logspace:
        b = (u1min-u1max)/np.log(u1min/u1max)
        a = u1max/b
        u1tab = [b*(a+np.log(x))
                 for x in np.linspace(u1min/u1max, 1,
                                      num_u1_steps-1)]
    else:
        u1tab = np.linspace(u1min, u1max, num_u1_steps-1).tolist()

    u1tab.append(u1ph)
    # Note: first num_u1_steps-1 are noload operation only
    parvardef = {
        "decision_vars": [
            {"values": u1tab, "name": "u1"},
            {"values": f1tab, "name": "f1"},
            {"values": [f/p for f in f1tab[:-1]] + [(1-slip)*f1tab[-1]/p],
             "name": "speed"}
        ]
    }

    parvar = femagtools.parstudy.List(
        workdir, condMat=condMat,
        magnetizingCurves=magnetizingCurves)

    # set AC simulation
    Q2 = machine['rotor']['num_slots']
    da1 = machine['bore_diam']
    slotmodel = [k for k in machine['rotor'] if isinstance(
        machine['rotor'][k], dict)][-1]
    Dr = (da1 - 2*machine['airgap'] -
          machine['rotor'][slotmodel].get('slot_height', 0) -
          machine['rotor'][slotmodel].get('slot_h1', 0))
    bar_len = machine['lfe']+np.pi*Dr/Q2/np.sin(np.pi*p/Q2)

    simulation = dict(
        calculationMode="asyn_motor",
        bar_len=bar_len,
        wind_temp=20,
        bar_temp=20,
        speed=f1/p,
        f1=f1,
        airgap_induc=True,
        num_par_wdgs=machine['windings'].get('num_par_wdgs', 1),
        wdgcon=CON[wdgcon],  # 0:open, 1:star, 2:delta
        u1=u1ph)  # phase voltage

    # start calculation
    m = copy.deepcopy(machine)
    try:
        m['windings'].pop('resistance')
    except KeyError:
        pass
    results = parvar(parvardef, m, simulation, engine)

    if simulation['wdgcon'] == 1:
        connu, conni = 1/np.sqrt(3), 1
    elif simulation['wdgcon'] == 2:
        connu, conni = 1, 1/np.sqrt(3)
    else:
        connu, conni = 1, 1
    u1tab = [connu*f['u1'][0] for f in results['f']]
    i1 = [conni*f['i1'][0] for f in results['f']]
    # amplitude of flux density in airgap
    bamp = [f['airgap']['Bamp'] for f in results['f']]
    wdg = femagtools.windings.Winding(
        {'Q': machine['stator']['num_slots'],
         'p': machine['poles']//2,
         'm': machine['windings']['num_phases'],
         'l': machine['windings']['num_layers'],
         'yd': machine['windings']['coil_span']})
    taup = np.pi*da1/2/p
    lfe = machine['lfe']
    n1 = wdg.turns_per_phase(machine['windings']['num_wires'],
                             machine['windings']['num_par_wdgs'])
    # main flux per phase at no load in airgap
    psih = [n1*wdg.kw()*taup*lfe*np.sqrt(2)/np.pi*b
            for b in bamp[:-1]]
    w1 = 2*np.pi*results['f'][0]['f1'][0]
    u1ref = u1ph
    psiref = float(ip.interp1d(u1tab[:-1], psih, kind='cubic')(u1ref))

    fitp, cov = so.curve_fit(
        lambda x, iml, ims, mexp: iml*x/psiref + ims*(x/psiref)**mexp,
        psih, i1[:-1], (1, 1, 1))
    iml, ims, mexp = fitp

    try:
        r1 = machine['windings']['resistance']
    except KeyError:
        n = machine['windings']['num_wires']
        g = simulation['num_par_wdgs']
        hs = machine['stator'].get('slot_height', 0)
        if 'dia_wire' in machine['windings']:
            aw = np.pi*machine['windings']['dia_wire']**2/4
        # TODO: read nc file and get slot area:
        # as = nc.windings[0].subregions[0].area()
        # aw = as/wdg.l/n*fillfac
        sigma = 56e6
        if 'material' in machine['winding']:
            sigma = condMat[machine['winding']]['elconduct']
        r1 = wdg_resistance(
            wdg, n, g, aw, da1, hs, lfe, sigma)

    # must correct psiref
    lh = psiref/(iml+ims)
    ls1 = results['f'][-1]['lh'] - (1+wdg.harmleakcoeff())*lh

    def imag(x):
        return iml*x/psiref + ims*(x/psiref)**mexp
    psix = so.fsolve(
        lambda psi: u1ref - np.sqrt(
            (r1*imag(psi))**2 +
            (w1*ls1 + w1*psi)**2),
        psiref)[0]
    logger.info(psix)
    psiref = psix
    fitp, cov = so.curve_fit(
        lambda x, iml, ims, mexp: iml*x/psiref + ims*(x/psiref)**mexp,
        psih, i1[:-1], (1, 1, 1))
    iml, ims, mexp = fitp

    p = results['f'][0]['p']
    m = results['f'][0]['num_phases']
    lh = psiref/(iml+ims)
    ls1 = results['f'][-1]['lh'] - lh
    r2 = results['f'][-1]['r2']
    ls2 = results['f'][-1]['ls2']
    pfe = results['f'][-1]['pfe1'][0]
    rotor_mass = sum([results['f'][-1].get('conweight', 0),
                      results['f'][-1].get('lamweight', 0)])

    n = machine['windings']['num_wires']
    g = simulation['num_par_wdgs']
    ag = machine['airgap']
    return {
        'p': p, 'm': m,
        'f1ref': f1, 'u1ref': u1ph,
        'rotor_mass': rotor_mass, 'kfric_b': 1,
        'r1': r1, 'r2': r2,
        'lsigma1': ls1, 'lsigma2': ls2,
        'psiref': psiref, 'wref': w1,
        'fec': pfe, 'fee': 0, 'fexp': 7.0,
        'iml': iml, 'ims': ims, 'mexp': mexp}


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
