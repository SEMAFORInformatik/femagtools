"""
  femagtools.im
  ~~~~~~~~~~~~~

  induction machine (electrical circuit model)

  Copyright 2022: Semafor Informatik & Energie AG, Switzerland
"""
import numpy as np
import scipy.interpolate as ip
import scipy.optimize as so
import logging
import pathlib
from .utils import skin_resistance, skin_leakage_inductance, wdg_leakage_inductances
import femagtools.windings
import femagtools.parstudy
import json
import copy
import warnings
import time

EPS = 1e-13


def log_interp1d(x, y, kind='cubic'):
    """logarithmic interpolation function, including extrapolation

    due to the nature of logarithmic functions, it is not possible to have
    a supporting point at 0!

    arguments:
    x, y (arraylike): supporting points
    kind (str): kind of underlying interpolation, default = 'cubic'
    """
    xx = np.asarray(x)
    yy = np.asarray(y)
    # TODO: check for x=0, dimensions, etc.
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = ip.interp1d(logx, logy, kind=kind, fill_value="extrapolate")
    def log_interp(zz): return np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp


eecdefaults = dict(
    zeta1=0.2,
    zeta2=2.4,
    gam=0,  # gam=0.7, 0..1
    kh=1,  # kh=2, number of vertical conductors in slot
    tcu1=20,
    rotor_mass=0,
    kfric_b=1,
    pl2v=0.5,
    tcu2=20)

logger = logging.getLogger('im')
logging.captureWarnings(True)


def ring_leakage_inductance(machine):
    """returns the ring leakage inductance
    ref: Design of Rotating Electrical Machines
    Juha Pyrhönen, Tapani Jokinen, Valeria Hrabovcova
    (Ed. 2008)  page 249
    """
    mue0 = 4*np.pi*1e-7

    Qr = machine['rotor']['num_slots']
    m = machine['windings']['num_phases']
    p = machine['poles']//2
    lbar = machine['lfe']
    ls = lbar
    nue = 0.36 if p == 1 else 0.18
    da1 = machine['bore_diam']
    ag = machine['airgap']
    slotmodel = [k for k in machine['rotor'] if isinstance(
        machine['rotor'][k], dict)][-1]

    Dr = (da1 - 2*ag -
          machine['rotor'][slotmodel].get('slot_height', 0) -
          machine['rotor'][slotmodel].get('slot_h1', 0))
    return mue0*Qr/m/p**2/3*((lbar-ls) + nue*np.pi*Dr/2/p)


def slot_opening_factor(n, p, bs, D, Q):
    k = 1 - bs/(D*np.pi/Q)
    return np.sin(n*p*np.pi/Q*(1-k))/(n*p*np.pi/Q*(1-k))


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
        self.kth1 = self.kth2 = 0.0039  # temp coefficient of stator wdg and rotor bar
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
            def imag(self, psi):
                """magnetizing current"""
                return psi/self.lh
            self._imag = imag
        elif 'iml' in parameters:
            def imag(psi):
                return (self.iml * np.abs(psi)/self.psiref +
                        self.ims*np.power(np.abs(psi)/self.psiref, self.mexp))
            self._imag = imag
        elif 'im' in parameters:
            self._imag = log_interp1d(self.psi, self.im)
            self.psi = self.psiref
        if 'pfe' in parameters:
            self.rh = self.m*(self.psiref*self.wref)**2/self.pfe
        for k in eecdefaults.keys():
            if not hasattr(self, k):
                setattr(self, k, eecdefaults[k])

    def imag(self, psi):
        """magnetizing current"""
        if np.isscalar(psi):
            return float(self._imag(psi))
        return self._imag(psi)

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

    def lrot(self, w):
        """rotor leakage inductance"""
        return skin_leakage_inductance(self.lsigma2, w, self.tcu2,
                                       self.zeta2, 1, self.pl2v)

    def lstat(self, w):
        """stator leakage inductance"""
        return self.lsigma1

    def rrot(self, w):
        """rotor resistance"""
        return skin_resistance(self.r2, w, self.tcu2, self.zeta2,
                               0.0, 1, kth=self.kth2)

    def rstat(self, w):
        """stator resistance"""
        return skin_resistance(self.r1, w, self.tcu1, self.zeta1,
                               self.gam, self.kh, kth=self.kth1)
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
        imag = complex(self.imag(psi))
        if abs(w1) > 0:
            imag += w1*psi/self.rfe(w1, psi)*1j
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
        # sk = self.rrot(0.)/(wsync*(self.lstat(wsync) + self.lrot(0.0)))
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

    def torque_chart(self, smin=-0.1, smax=0.1, nsamples=100):
        """
        calculate torque(s) curve

        :param smin: min slip
        :param smax: max slip
        :return: dict with slip and torque lists
        """

        s = np.linspace(smin, smax, nsamples)
        f1 = self.f1ref
        u1 = self.u1ref
        w1 = 2 * np.pi * f1
        wm = (1-s) * w1 / self.p
        r = {}
        r['s'] = list(s)
        r['T'] = [self.torqueu(w1, u1, wx) for wx in wm]
        return r

    def operating_point(self, T, n, u1max, Tfric=None):
        """
        calculate single operating point.

        return dict with values for
            s -- (float) slip at this op
            f1 -- (float) stator frequency in Hz
            u1 -- (float) stator phase voltage in V rms
            i1 -- (float) stator phase current in A rms
            i2 -- (float) rotor bar current in A rms
            p1 -- (float) electrical power in W
            pmech -- (float) mechanical power in W
            plcu1 -- (float) stator copper losses in W
            plcu2 -- (float) rotor copper losses in W
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
            Tfric -- (float, optional) the friction torque to consider in Nm.
                     If Tfric is None, the friction torque is calculated according to the rotor-mass and the
                     frition factor of the machine (kfric_b). Friction is written to the result dict!
        """

        r = {}  # result dit

        if Tfric:
            tfric = Tfric
        else:
            # formula from
            # PERMANENT MAGNET MOTOR TECHNOLOGY: DESIGN AND APPLICATIONS
            # Jacek Gieras
            #
            # plfric = kfric_b m_r n 1e-3 W
            # -- kfric_b : 1..3 W/kg/rpm
            # -- m_r: rotor mass in kg
            # -- n: rotor speed in rpm
            tfric = self.kfric_b * self.rotor_mass * 30e-3 / np.pi
            # TODO: make frictiontorque speed depended?

        tq = T + tfric
        wm = 2*np.pi*n
        w1 = self.w1(u1max, self.psiref, tq, wm)
        s = (w1 - self.p*wm) / w1
        r['s'] = float(s)
        r['f1'] = float(w1/2/np.pi)
        u1 = self.u1(w1, self.psi, wm)
        r['u1'] = float(np.abs(u1))
        i1 = self.i1(w1, self.psi, wm)
        r['i1'] = float(np.abs(i1))
        r['cosphi'] = float(np.cos(np.angle(u1) - np.angle(i1)))
        r['plfe1'] = float(self.m * np.abs(u1) ** 2 / self.rfe(w1, self.psi))
        i2 = self.i2(w1, self.psi, wm)
        r['i2'] = float(np.abs(i2))
        r['plcu1'] = float(self.m * np.abs(i1) ** 2 * self.rstat(w1))
        r['plcu2'] = float(self.m * np.abs(i2) ** 2 *
                           self.rrot(w1 - self.p * wm))
        r['plfric'] = float(2 * np.pi * n * tfric)
        pmech = 2 * np.pi * n * tq
        r['pmech'] = float(pmech)
        pltotal = r['plfe1'] + r['plfric'] + r['plcu1'] + r['plcu2']
        r['losses'] = float(pltotal)
        p1 = pmech + pltotal
        r['p1'] = float(p1)
        if np.abs(pmech) < 1e-12:
            eta = 0  # power to low for eta calculation
        elif p1 > pmech:
            eta = pmech/p1  # motor
        else:
            eta = p1/pmech  # generator
        r['eta'] = float(eta)
        r['T'] = float(T)
        r['Tfric'] = float(tfric)
        r['n'] = float(n)
        r['tcu1'] = float(self.tcu1)
        r['tcu2'] = float(self.tcu2)

        return r

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

        r = dict(u1=[], i1=[], T=[], cosphi=[], n=[], s=[], sk=[],
                 plfe1=[], plcu1=[], plcu2=[])
        T = [tload2(wx) for wx in wmtab]
        tfric = self.kfric_b*self.rotor_mass*30e-3/np.pi
        w1tab = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
                r['s'].append(float((w1 - self.p * wm) / w1))
                r['sk'].append(self.sk(w1, u1/w1))
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
    f1 -- line voltage frequency [Hz]
    u1 -- line voltage (line-to-line) [V]
    wdgcon -- winding connection  ('open', 'wye', 'star', 'delta')
    machine -- dict() with machine parameters
    magnetizingCurves -- list of dict() with BH curves (or directory with MC/MCV files)
    condMat -- list of dict() with conductor material properties

    optional:
    num_steps: (int) number of current steps for the no-load
        calculation (default 5)
    logspace: (bool) uses log distributed current samples if true
         (linear otherwise) (default true)
    i_max_fact: (float) factor for maximum current to calculate no_load flux (default=2.5)
    templatedirs: (list of str) names of directories to search for templates
    """
    CON = {'open': 0, 'wye': 1, 'star': 1, 'delta': 2}
    p = machine['poles']//2
    slip = 1e-2
    u1ph = u1
    w1 = 2*np.pi*f1
    if CON[wdgcon] == 1:
        u1ph /= np.sqrt(3)
    wdg = femagtools.windings.Winding(
        {'Q': machine['stator']['num_slots'],
         'm': machine['windings']['num_phases'],
         'p': machine['poles']//2,
         'l': machine['windings']['num_layers'],
         'yd': machine['windings']['coil_span']})
    L1 = wdg.inductance(
        nwires=machine['windings']['num_wires'],
        g=machine['windings']['num_par_wdgs'],
        da1=machine['bore_diam'],
        lfe=machine['lfe'],
        ag=machine['airgap'])
    i1max_fact = kwargs.get('i_max_fact', 2.5)
    i1max = i1max_fact*u1ph/w1/L1
    i1min = u1ph/w1/L1/5
    num_steps = kwargs.get('num_steps', 5)
    if kwargs.get('logspace', False):
        b = (i1min-i1max)/np.log(i1min/i1max)
        a = i1max/b
        i1tab = [b*(a+np.log(x))
                 for x in np.linspace(i1min/i1max, 1,
                                      num_steps)]
    else:
        i1tab = np.linspace(i1min, i1max, num_steps).tolist()

    m = copy.deepcopy(machine)
    Q2 = m['rotor']['num_slots']
    noloadsim = dict(
        calculationMode="noloadflux-rot",
        curvec=i1tab,
        num_par_wdgs=machine['windings'].get('num_par_wdgs', 1),
        Q2=Q2)

    da1 = m['bore_diam']
    ag = m['airgap']
    # do not create airgap nodechains automatically
    # we prefer an integer number of elements per rotor slot
    nper = np.lcm(wdg.Q, Q2)/Q2
    m['num_agnodes'] = Q2*np.ceil(2*(da1-ag)*np.pi/Q2/ag/nper)*nper
    try:
        m['windings'].pop('resistance')
    except KeyError:
        pass

    parstudy = femagtools.parstudy.ParameterStudy(
        workdir, condMat=condMat,
        magnetizingCurves=magnetizingCurves)

    builder = femagtools.fsl.Builder(kwargs.get('templatedirs', []))
    model = femagtools.model.MachineModel(m)
    # modelfiles = parstudy.setup_model(builder, model)
    extra_result_files = [f'noloadbag-{i+1}.dat'
                          for i in range(num_steps)] + [
        'psi-rot-mag.dat']

    # set AC simulation
    slotmodel = [k for k in machine['rotor'] if isinstance(
        machine['rotor'][k], dict)][-1]
    Dr = (da1 - 2*ag -
          machine['rotor'][slotmodel].get('slot_height', 0) -
          machine['rotor'][slotmodel].get('slot_h1', 0))
    bar_len = machine['lfe']+np.pi*Dr/Q2/np.sin(np.pi*p/Q2)

    ecsim = dict(
        calculationMode="ec-rotorbar",
        bar_len=bar_len,
        wind_temp=20,
        bar_temp=20,
        f1=f1)

    rotorbar = copy.deepcopy(machine)
    # remove stator slot model and windings
    for k in rotorbar['stator']:
        if isinstance(rotorbar['stator'][k], dict):
            d = rotorbar['stator'].pop(k)
            break
    d = rotorbar.pop('windings')
    if 'shaft_diam' in rotorbar:
        rotorbar['inner_diam'] = rotorbar.pop('shaft_diam')
    # use one rotor slot only
    rotorbar['stator']['statorRing'] = {'num_slots': Q2}
    for k in rotorbar['rotor']:
        if isinstance(rotorbar['rotor'][k], dict):
            rotorbar['rotor'][k]['num_slots_gen'] = 1
            rotorbar['rotor'][k]['zeroangle'] = 90-180/Q2
            break

    loadsim = dict(  # not used
        calculationMode="asyn_motor",
        bar_len=bar_len,
        wind_temp=20,
        bar_temp=20,
        speed=(1-slip)*f1/p,
        f1=f1,
        num_par_wdgs=machine['windings'].get('num_par_wdgs', 1),
        wdgcon=CON[wdgcon],  # 0:open, 1:star, 2:delta
        u1=u1ph)  # phase voltage

    # prepare calculation
    job = engine.create_job(workdir)
    task = job.add_task(_eval_noloadrot(), extra_result_files)
    logger.debug("Task %s noload workdir %s result files %s",
                 task.id, task.directory, task.extra_result_files)
    # create model
    for mc in parstudy.femag.copy_magnetizing_curves(
            model,
            dir=task.directory):
        task.add_file(mc)

    task.add_file(
        'femag.fsl',
        builder.create_model(model, condMat=parstudy.femag.condMat) +
        builder.create_analysis(noloadsim) +
        ['save_model("close")'])

    # ec simulation
    barmodel = femagtools.model.MachineModel(rotorbar)
    extra_result_files = ['bar.dat']
    r = (da1-ag)/2
    task = job.add_task(_eval_ecsim())
    logger.debug("Task %s rotobar workdir %s result files %s",
                 task.id, task.directory, task.extra_result_files)
    task.set_stateofproblem('mag_dynamic')
    for mc in parstudy.femag.copy_magnetizing_curves(
            barmodel,
            dir=task.directory,
            recsin='cur'):  # todo cur
        task.add_file(mc)
    task.add_file(
        'femag.fsl',
        builder.create_model(barmodel,
                             condMat=parstudy.femag.condMat) +
        builder.create_analysis(ecsim) +
        ['save_model("close")'])

    # AC simulation
    actask = job.add_task(result_files=['end_wind_leak.dat'])
    logger.debug("Task %s loadsim workdir %s result files %s",
                 task.id, task.directory, task.extra_result_files)
    actask.set_stateofproblem('mag_dynamic')
    for mc in parstudy.femag.copy_magnetizing_curves(
            model,
            dir=actask.directory,
            recsin='flux'):
        actask.add_file(mc)
    actask.add_file(
        'femag.fsl',
        builder.create_model(model, condMat=parstudy.femag.condMat) +
        builder.create_analysis(loadsim) +
        ['save_model("close")'])
    # start FE simulations
    tstart = time.time()
    status = engine.submit()
    logger.info('Started %s', status)
    status = engine.join()
    tend = time.time()
    logger.info("Elapsed time %d s Status %s",
                (tend-tstart), status)
    # collect results
    results = []
    errors = []
    for t in job.tasks:
        if t.status == 'C':
            results.append(t.get_results())
        else:
            logger.warning("Status %s", t.status)
            results.append({})

    i1_0 = results[0]['i1_0'].tolist()
    psi1_0 = results[0]['psi1_0'].tolist()
    # amplitude of flux density in airgap
    bamp = results[0]['Bamp']
    taup = np.pi*(da1-ag)/2/p
    lfe = machine['lfe']
    n1 = wdg.turns_per_phase(machine['windings']['num_wires'],
                             machine['windings']['num_par_wdgs'])
    # main flux per phase at no load in airgap
    kfe = machine['stator'].get('fillfac', 1.0)
    psih = np.mean([[2/np.pi*n1*wdg.kw()*taup*lfe*b/np.sqrt(2)
                   for b in bb]
                    for bb in bamp], axis=1)
    i10tab = [0] + i1_0
    psihtab = [0] + psih.tolist()
    u1ref = u1ph
    psiref = u1ref/w1

    # def inoload(x, iml, ims, mexp):
    #    """return noload current"""
    #    return iml*x/psiref + ims*(x/psiref)**mexp
    # fitp, cov = so.curve_fit(inoload, psihtab, i10tab, (1, 1, 1))
    # iml, ims, mexp = fitp
    # logger.info("iml, ims, mexp %g, %g, %g",
    #            iml, ims, mexp)
    # i1tab.insert(0, 0)
    # psi1_0.insert(0, 0)
    # i1_0.insert(0, 0)
    logger.info("psi1_0 %s", np.mean(psi1_0, axis=1))
    logger.info("psih %s", psih)
    logger.debug("psi1_0-psih %s", np.mean(psi1_0, axis=1)-psih)
    logger.debug("i1tab %s", i1tab)
    logger.debug("u1ref %s", u1ref)
    logger.debug("w1 %s", w1)
    logger.debug("L1 %s", L1)
    try:
        r1 = machine['windings']['resistance']
    except KeyError:
        from .utils import wdg_resistance
        slotmodel = [k for k in machine['stator'] if isinstance(
            machine['stator'][k], dict)][-1]
        if slotmodel == 'stator1':
            hs = machine['stator']['stator1']['slot_rf1'] - \
                machine['stator']['stator1']['tip_rh1']
        else:
            hs = machine['stator'][slotmodel].get('slot_height',
                                                  0.33*(machine['outer_diam']-da1))
        n = machine['windings']['num_wires']
        if 'dia_wire' in machine['windings']:
            aw = np.pi*machine['windings'].get('dia_wire', 1e-3)**2/4
        else:  # wire diameter from slot area
            aw = 0.75 * \
                machine['windings'].get(
                    'cufilfact', 0.45)*np.pi*da1*hs/wdg.Q/2/n
        # TODO: read nc file and get slot area:
        # as = nc.windings[0].subregions[0].area()
        # aw = as/wdg.l/n*fillfac
        # TODO: sigma = 58e6
        # if 'material' in machine['windings']:
        #    sigma = condMat[machine['windings']]['elconduct']
        g = loadsim['num_par_wdgs']
        r1 = wdg_resistance(
            wdg, n, g, aw, da1, hs, lfe)

    # psi1 = ip.interp1d(i1_0, np.mean(psi1_0, axis=1),
    #                    kind='quadratic')
    psi1 = log_interp1d(i1_0, np.mean(psi1_0, axis=1))
    imref = so.fsolve(
        lambda imx: u1ref - np.sqrt((imx*r1)**2 +
                                    (w1*psi1(imx))**2),
        u1ph/w1/L1)[0]
    # psihref = float(ip.interp1d(i1tab, psih,
    #                             kind='quadratic')(imref))
    psihref = float(log_interp1d(i1tab, psih)(imref))
    psi1ref = float(psi1(imref))

    lh = psihref/imref
    L1 = psi1ref/imref
    ls1 = L1 - lh*(1+wdg.harmleakcoeff())
    logger.info('L1: %g Lh: %g Ls1: %g kw %g, sigma0: %g',
                L1, lh, ls1, wdg.kw(), wdg.harmleakcoeff())
    ü = 3*4*(wdg.kw()*n1)**2/Q2
    r2 = results[1]['r2']*ü
    ls2 = results[1]['ls2']*ü + ring_leakage_inductance(machine)*ü
    zeta2 = results[1]['zeta2']
    pl2v = results[1]['pl2v']

    # end winding leakage
    basedir = pathlib.Path(actask.directory)
    leakfile = basedir / 'end_wind_leak.dat'
    try:
        leakages = [float(x)
                    for x in leakfile.read_text().split()]
        ls1 += leakages[1]  # TODO: np.linalg.norm(leakages[1:])
        logger.info("L1ew %s", leakages[1])
    except:
        logger.warning("No end winding leakage")
        Lu, Lew = wdg_leakage_inductances(machine)
        logger.info("L1ew %s", Lew)
        ls1 += Lew
    p = results[2]['p']
    m = results[2]['num_phases']
    # r2 = results[1]['r2']
    # ls2 = results[1]['ls2']
    pfe = results[2]['pfe1'][0]
    rotor_mass = sum([results[2].get('conweight', 0),
                      results[2].get('lamweight', 0)])

    n = machine['windings']['num_wires']
    g = loadsim['num_par_wdgs']
    impars = {
        'p': p, 'm': m,
        'f1ref': f1, 'u1ref': u1ph,
        'rotor_mass': rotor_mass, 'kfric_b': 1,
        'r1': r1, 'r2': r2,
        'lsigma1': ls1,
        'lsigma2': ls2, 'zeta2': zeta2, 'pl2v': pl2v,
        'psiref': psihref, 'wref': w1,
        'fec': pfe, 'fee': 0, 'fexp': 7.0,
        'im': i1tab, 'psi': psih.tolist()}
    try:
        wmat = machine['windings']['material']
        impars['kth1'] = parstudy.femag.condMat.find(wmat)['tempcoef']
    except KeyError:
        logger.warning('Missing winding material id')
    try:
        bmat = machine['rotor']['material']
        impars['kth2'] = parstudy.femag.condMat.find(bmat)['tempcoef']
    except KeyError:
        logger.warning('Missing winding material id')

    return impars


class _eval_noloaddc():
    """ Result Functor for noloadflux dc calc"""

    def __init__(self):
        pass

    def __call__(self, task):
        basedir = pathlib.Path(task.directory)
        psimag = np.loadtxt(basedir/'noloadflux.dat').T
        im, psi = psimag[0:6:2], psimag[12::2]
        d = 0
        a = np.array(
            (np.cos(d), np.cos(d-2*np.pi/3), np.cos(d+2*np.pi/3)))/3*2
        i0 = np.sum(a*im.T, axis=1)/np.sqrt(2)
        psi0 = np.sum(a*psi.T, axis=1)/np.sqrt(2)

        bag = np.loadtxt(basedir / 'noloadbag.dat').T
        Bamp = [b['Bamp'] for b in [
            femagtools.airgap.fft(bag[0], b)
            for b in bag[1:]]]
        return {'i1_0': i0, 'psi1_0': psi0, 'Bamp': Bamp}


class _eval_noloadrot():
    """ Result Functor for noloadflux rot calc"""

    def __init__(self):
        pass

    def __call__(self, task):
        basedir = pathlib.Path(task.directory)
        psimag = np.loadtxt(basedir/'psi-rot-mag.dat')
        pos = np.unique(psimag[:, 0])
        ncurs = psimag[:, 0].shape[0]//pos.shape[0]
        ire = psimag[:ncurs, 1:7:2]  # assuming same current for all steps
        psire = np.reshape(psimag[:, 7:12:2], (-1, ncurs, 3))

        i0 = np.linalg.norm(
            femagtools.machine.T(0).dot(ire.T),
            axis=0)/np.sqrt(2)

        psi0 = np.array([np.linalg.norm(
            femagtools.machine.T(0).dot(psire[:, k, :].T),
            axis=0)/np.sqrt(2)
            for k in range(ncurs)])

        # matrix (i x j x k) of curr, rotor pos, angle
        Bamp = [[femagtools.airgap.fft(bags[:, 0], b)['Bamp']
                 for b in bags.T[1:]]
                for bags in [np.loadtxt(p) for p in sorted(
                    basedir.glob(f"noloadbag-*.dat"),
                    key=lambda path: int(path.stem.rsplit("-", 1)[1]))]]
        return {'i1_0': i0, 'psi1_0': psi0, 'Bamp': Bamp}


class _eval_ecsim():
    """ Result Functor for ec simulation"""

    def __init__(self):
        pass

    def __call__(self, task):
        import lmfit
        from .utils import xiskin, kskinl
        basedir = pathlib.Path(task.directory)
        freqrind = np.loadtxt(basedir/'bar.dat').T
        rbar = freqrind[1]
        lbar = freqrind[2:]
        f = freqrind[0]
        r0 = rbar[0]
        l0 = lbar[0, 0]
        temp = 20

        def barimp(f, zeta, pl2v):
            w = 2*np.pi*f
            r = skin_resistance(r0, w, temp, zeta)
            xi = xiskin(w, temp, zeta)
            lbar = l0*(1 + pl2v*(kskinl(xi, 1)-1))
            return r + 1j*lbar
        model = lmfit.model.Model(barimp)
        params = model.make_params(zeta=1, pl2v=0.5)
        guess = lmfit.models.update_param_vals(params, model.prefix)
        imp = rbar + 1j*lbar[-1]
        rfit = model.fit(imp, params=guess, f=f, verbose=True)
        zeta = rfit.params['zeta'].value
        pl2v = rfit.params['pl2v'].value
        logger.info("ECSIM %s",
                    {'r2': r0, 'ls2': l0, 'zeta2': zeta, 'pl2v': pl2v})
        return {'r2': r0, 'ls2': l0, 'zeta2': zeta, 'pl2v': pl2v}


class _eval_end_wind_leak():
    """ Result Functor for ac simulation"""

    def __init__(self):
        pass

    def __call__(self, task):
        basedir = pathlib.Path(task.directory)
        leakfile = basedir / 'end_wind_leak.dat'
        if leakfile.exists():
            leakages = leakfile.read_text().split()
            return {'lse': np.linalg.norm(leakages[1:])}
        return {'lse': 0}


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
