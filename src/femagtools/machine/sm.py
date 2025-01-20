"""Wound-rotor synchronous machine (EESM) electrical circuit model

"""
import logging
import warnings
import pathlib
import numpy as np
import scipy.optimize as so
import scipy.interpolate as ip
from .utils import skin_resistance, wdg_resistance, betai1, iqd, KTH, create_wdg
from .. import parstudy, windings
import femagtools.bch

EPS = 1e-13

eecdefaults = {
    'zeta1': 0.3,
    'zeta2': 0,
    'gam': 0.7,
    'kh': 2,
    'tcu1': 20,
    'tcu2': 20,
    'rotor_mass': 0,
    'kfric_b': 1,
    'kpfe': 1 # iron loss factor
}

logger = logging.getLogger('sm')
logging.captureWarnings(True)


def parident(workdir, engine, machine,
             magnetizingCurves, condMat,
             **kwargs):
    """return dict of parameters of equivalent circuit for
        electrically excited synchronous machines

        arguments:
        machine: dict() with machine parameters
        magnetizingCurves: list of dict() with BH curves
        condMat: list of dict() with conductor material properties

        optional arguments:
        num_cur_steps: number of current steps (default 5)
        num_beta_steps: number of current steps (default 13)
        num_exc_steps: number of excitation current (default 7)
        speed: rotor speed in 1/s (default 160/p)
        i1_max: maximum current in A rms (default approx 3*i1nom)
        beta_min: minimal current angle (default -180°)
        beta_max: maximal current angle (default 0°)
        num_move_steps: number of move steps
        cmd: femag command (default None, platform executable)
        """
    cmd = kwargs.get('cmd', None)

    wdgk = 'windings' if 'windings' in machine else 'winding'
    g = machine[wdgk].get('num_par_wdgs', 1)
    N = machine[wdgk]['num_wires']
    if 'cufilfact' in machine[wdgk]:
        fcu = machine[wdgk]['cufilfact']
    elif 'fillfac' in machine[wdgk]:
        fcu = machine[wdgk]['fillfac']
    else:
        fcu = 0.42
    try: # calc basic dimensions if not fsl or dxf model
        from ..model import MachineModel
        wdg = create_wdg(machine)
        Q1 = wdg.Q
        model = MachineModel(machine)
        Jmax = 20e6  # max current density in A/m2
        Acu = fcu*model.slot_area()  # approx. copper area of one slot
        i1_max = round(g*Acu/wdg.l/N*Jmax/10)*10
    except KeyError:
        if kwargs.get('i1_max', 0) == 0:
            raise ValueError('i1_max missing')
        i1_max = kwargs['i1_max']

    ifnom = machine['rotor']['ifnom']
    exc_logspace = True
    ifmin, ifmax = ifnom/4, 1.4*ifnom
    if exc_logspace:
        excur = np.logspace(np.log(ifmin), np.log(ifmax),
                            kwargs.get("num_exc_steps", 6),
                            base=np.exp(1)).tolist()
    else:
        excur = np.linspace(ifmin, ifmax,
                            kwargs.get("num_exc_steps", 6))

    logger.info("Exc current %s", excur)
    parvardef = {
        "decision_vars": [
            {"values": excur, "name": "load_ex_cur"}
        ]
    }

    parvar = parstudy.List(
        workdir,  condMat=condMat,
        magnetizingCurves=magnetizingCurves, cmd=cmd)

    simulation = dict(
        calculationMode=kwargs.get('calculationMode',
                                   'ld_lq_fast'),
        wind_temp=20.0,
        i1_max=kwargs.get('i1_max', i1_max),
        maxid=kwargs.get('maxid', 0),
        minid=kwargs.get('minid', -i1_max),
        maxiq=kwargs.get('maxiq', i1_max),
        miniq=kwargs.get('miniq', -i1_max),
        delta_id=kwargs.get('delta_id', i1_max/5),
        delta_iq=kwargs.get('delta_iq', i1_max/5),
        beta_min=kwargs.get('beta_min', -180),
        beta_max=kwargs.get('beta_max', 0),
        num_move_steps=kwargs.get('num_move_steps', 31),
        load_ex_cur=0.5,
        num_cur_steps=kwargs.get('num_cur_steps', 5),
        num_beta_steps=kwargs.get('num_beta_steps', 13),
        num_par_wdgs=machine[wdgk].get('num_par_wdgs', 1),
        skew_angle=kwargs.get('skew_angle', 0.0),
        num_skew_steps=kwargs.get('num_skew_steps', 0.0),
        period_frac=kwargs.get('period_frac', 6),
        speed=kwargs.get('speed', 50))

    ###self.cleanup()  # remove previously created files in workdir
    results = parvar(parvardef, machine, simulation, engine)
    b = results['f'][-1]

    if 'poles' not in machine:  # dxf model?
        machine['poles'] = 2*results['f'][0]['machine']['p']
        da1 = 2*results['f'][0]['machine']['fc_radius']
        wdg = create_wdg(machine)

    if simulation['calculationMode'] == 'ld_lq_fast':
        idname = 'ldq'
    else:
        idname = 'psidq'
    for r in results['f']:
        for k in ('hf', 'ef'):
            if k in r['lossPar']:
                r[idname]['losses'][k] = r['lossPar'][k]
    try:
        rotor_mass = sum(results['f'][-1]['weights'][-1])
    except KeyError:
        rotor_mass = 0  # need femag classic > rel-9.3.x-48-gca42bbd0

    losskeys = ('speed', 'ef', 'hf', 'cf',
                'styoke_hyst', 'stteeth_hyst', 'styoke_eddy',
                'stteeth_eddy', 'rotor_hyst', 'rotor_eddy',
                'styoke_excess', 'stteeth_excess', 'rotor_excess')

    # winding resistance
    try:
        r1 = machine[wdgk]['resistance']
    except KeyError:
        try:
            Q1 = machine['stator']['num_slots']
            yd = machine[wdgk].get('coil_span', Q1/machine['poles'])
            wdg = windings.Winding(
                {'Q': Q1,
                 'm': machine[wdgk]['num_phases'],
                 'p': machine['poles']//2,
                 'l': machine[wdgk]['num_layers'],
                 'yd': yd})

            lfe = machine['lfe']
            da1 = machine['outer_diam']
            if 'dia_wire' in machine[wdgk]:
                aw = np.pi*machine[wdgk].get('dia_wire', 1e-3)**2/4
            else:  # wire diameter from slot area
                aw = 0.75 * fcu * np.pi*da1*hs/Q1/wdg.l/N
            r1 = wdg_resistance(wdg, N, g, aw, da1, hs, lfe)
        except (KeyError, NameError):
            from .. import nc
            model = nc.read(str(pathlib.Path(workdir) / machine['name']))
            try:
                nlayers = wdg.l
            except UnboundLocalError:
                wdg = create_wdg(machine)
                nlayers = wdg.l
                da1 = 2*results['f'][0]['machine']['fc_radius']
            Q1 = wdg.Q
            istat = 0 if model.get_areas()[0]['slots'] else 1
            asl = model.get_areas()[istat]['slots']
            # diameter of wires
            aw = fcu*asl/Q1/nlayers/N
            hs = asl/(np.pi*da1/3)
            r1 = wdg_resistance(wdg, N, g, aw, da1, hs, lfe)

    if simulation['calculationMode'] == 'ld_lq_fast':
        dqpars = dict(m=3, p=b['machine']['p'],
                      r1=float(r1),
                      r2=machine['rotor'].get('resistance', 1),
                      rotor_mass=rotor_mass, kfric_b=1,
                      ldq=[dict(
                          ex_current=b['machine']['ex_current'],
                          i1=b['ldq']['i1'],
                          beta=b['ldq']['beta'],
                          psid=b['ldq']['psid'],
                          psiq=b['ldq']['psiq'],
                          torque=b['ldq']['torque'],
                          ld=b['ldq']['ld'],
                          lq=b['ldq']['lq'],
                          losses={k: b['ldq']['losses'][k]
                                    for k in losskeys if k in b['ldq']['losses']})
                         for b in results['f']])
    else:
        dqpars = dict(m=3, p=b['machine']['p'],
                      r1=r1,
                      r2=machine['rotor'].get('resistance', 1),
                      rotor_mass=rotor_mass, kfric_b=1,
                      psidq=[dict(
                          ex_current=b['machine']['ex_current'],
                          iq=b['psidq']['iq'],
                          id=b['psidq']['id'],
                          psidq_ldq=b['psidq_ldq'],
                          psid=b['psidq']['psid'],
                          psiq=b['psidq']['psiq'],
                          torque=b['psidq']['torque'],
                          losses={k: b['psidq']['losses'][k]
                                for k in losskeys if k in b['psidq']['losses']})
                       for b in results['f']])

    if 'current_angles' in results['f'][0]:
        dqpars['current_angles'] = results['f'][0]['current_angles']
    return dqpars

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
    d = np.diff(exc)
    m = np.mean(d)
    return np.max(np.abs(d**2 - m**2)) < 1e-9


def _gradient_respecting_bounds(bounds, fun, eps=1e-8):
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
    """ represent Synchronous machine with wound rotor (EESM)

    Arguments:
        eecpars: dict() electrical circuit parameters
        """
    def __init__(self, eecpars, **kwargs):
        self.kth1 = KTH
        self.kth2 = KTH
        self.skin_resistance = [None, None]
        self.kpmag = 1
        self.bertotti = False
        # here you can set user defined functions for calculating the skin-resistance,
        # according to the current frequency w. First function in list is for stator, second for rotor.
        # If None, the femagtools intern default implementation is used.
        # User defined functions need to have the following arguments:
        # - r0: (float) dc-resistance at 20°C
        # - w: (float)  current frequency in rad (2*pi*f)
        # - tcu: (float) conductor temperature in deg Celsius
        # - kth: (float) temperature coefficient (Default = 0.0039, Cu)
        for k in eecdefaults.keys():
            setattr(self, k, eecdefaults[k])

        for k in eecpars:
            if k not in ('ldq', 'psidq'):
                setattr(self, k, eecpars[k])

        for k in kwargs:
            setattr(self, k, kwargs[k])

        try:
            self.tfric = self.kfric_b*self.rotor_mass*30e-3/np.pi
        except AttributeError:
            self.tfric = 0

        self.fo = 50
        self.plexp = {'styoke_hyst': [1.0, 1.0],
                      'stteeth_hyst': [1.0, 1.0],
                      'styoke_eddy': [2.0,2.0],
                      'stteeth_eddy': [2.0,2.0],
                      'rotor_hyst': [1.0,1.0],
                      'rotor_eddy': [2.0, 2.0]}

    def _set_losspar(self, speed, ef, hf):
        self.fo = speed*self.p
        self.plexp = {'styoke_hyst': hf,
                      'stteeth_hyst': hf,
                      'styoke_eddy': ef,
                      'stteeth_eddy': ef,
                      'rotor_hyst': hf,
                      'rotor_eddy': ef}

        if self.bertotti:
            self.plexp.update({
                'styoke_excess': 1.5,
                'stteeth_excess':1.5,
                'rotor_excess': 1.5})


    def pfric(self, n):
        """friction and windage losses"""
        return 2*np.pi*n*self.tfric

    def rstat(self, w):
        """stator resistance"""
        if isinstance(self.zeta1, list):
            # polyfit from ac loss calculation
            freq = w/2/np.pi
            kr = self.zeta1[0]*freq**3 + self.zeta1[1]*freq**2 + \
                self.zeta1[2]*freq + self.zeta1[3]
            if isinstance(kr, list):
                kr = np.array(kr)
                kr[kr<1.0] = 1.0
            elif isinstance(kr, np.ndarray):
                kr[kr<1.0] = 1.0
            else:
                if kr < 1.0:
                    kr = 1.0
            return self.r1*(1 + self.kth1*(self.tcu1 - 20))*kr  # ref 20°C
        sr = self.skin_resistance[0]
        if sr is not None:
            return sr(self.r1, w, self.tcu1, kth=self.kth1)
        else:
            return skin_resistance(self.r1, w, self.tcu1, self.zeta1,
                                   self.gam, self.kh, kth=self.kth1)

    def rrot(self, w):
        """rotor resistance"""
        sr = self.skin_resistance[1]
        if sr is not None:
            return sr(self.r2, w, self.tcu2, kth=self.kth2)
        else:
            return skin_resistance(self.r2, w, self.tcu2, self.zeta2,
                                   0.0, 1, kth=self.kth2)

    def torque_iqd(self, iq, id, iex):
        "torque at q-d-current"
        psid, psiq = self.psi(iq, id, iex)
        return self.m*self.p/2*(psid*iq - psiq*id)

    def tloss_iqd(self, iq, id, iex, n):
        """return loss torque of d-q current, iron loss correction factor
            and friction windage losses"""
        if n > 1e-3:
            f1 = self.p*n
            plfe = self.kpfe * (self.iqd_plfe1(iq, id, iex, f1)
                                + self.iqd_plfe2(iq, id, iex, f1))
            return (plfe + self.pfric(n))/(2*np.pi*n)
        return 0

    def tmech_iqd(self, iq, id, iex, n):
        """return shaft torque of d-q current and speed"""
        return self.torque_iqd(iq, id, iex) - self.tloss_iqd(iq, id, iex, n)

    def torquemax(self, i1, iex):
        "returns maximum torque of i1, iex (nan if i1 out of range)"
        def torquei1b(b):
            return self.torque_iqd(*iqd(b[0], i1), iex)
        res = so.minimize(torquei1b, (0,))
        return res.fun

    def torquemin(self, i1, iex):
        "returns minimum torque of i1, iex (nan if i1 out of range)"
        def torquei1b(b):
            return self.torque_iqd(*iqd(b[0], i1), iex)
        res = so.minimize(torquei1b, (-np.pi/2,))
        return -res.fun

    def uqd(self, w1, iq, id, iex):
        """return uq, ud of frequency w1 and d-q current"""
        psid, psiq = self.psi(iq, id, iex)
        r1 = self.rstat(w1)
        return r1*iq + w1*psid, r1*id - w1*psiq

    def plcu1(self, iqde, w1):
        r1 = self.rstat(w1)
        return 3/2*r1*(iqde[0]**2 + iqde[1]**2)

    def plcu2(self, iqde, w1=0):
        r2 = self.rrot(0)
        return r2*iqde[2]**2

    def culoss(self, iqde, w1=0):
        return self.plcu1(iqde, w1) + self.plcu2(iqde, 0)

    def iqd_plcu1(self, iq, id, f1):
        return self.plcu1((iq, id), 2*np.pi*f1)

    def iqd_plcu2(self, iq, id, iex):
        return self.plcu2((iq, id, iex))

    def iqd_plfe2_(self, iq, id, f1):
        return np.zeros(np.asarray(iq).shape)

    def iqd_plmag(self, iq, id, f1):
        return np.zeros(np.asarray(iq).shape)

    def iqd_tmech(self, torque, n, disp=False, maxiter=500):
        """return currents for shaft torque with minimal losses"""
        if torque > 0:
            startvals = self.bounds[0][1]/2, 0, self.bounds[-1][1]
        else:
            startvals = -self.bounds[0][1]/2, 0, self.bounds[-1][1]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            def sqrtculoss(iqde):
                pcu = self.culoss(iqde)
                return pcu

            res = so.minimize(
                self.culoss, startvals, method='SLSQP',  # trust-constr
                bounds=self.bounds,
                #            jac=gradient_respecting_bounds(self.bounds, self.culoss),
                constraints=[
                    {'type': 'eq',
                     'fun': lambda iqd: self.tmech_iqd(*iqd, n) - torque}])
            #options={'disp': disp, 'maxiter': maxiter})
            if res['success']:
                return res.x

        logger.warning("%s: torque=%f %f, io=%s",
                       res['message'], torque, self.tmech_iqd(*startvals, n),
                       startvals)
        raise ValueError(res['message'])

    def iqd_torque(self, torque, disp=False, maxiter=500):
        """return currents for torque with minimal losses"""
        if torque > 0:
            startvals = self.bounds[0][1]/2, 0, self.bounds[-1][1]
        else:
            startvals = -self.bounds[0][1]/2, 0, self.bounds[-1][1]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            def sqrtculoss(iqde):
                pcu = self.culoss(iqde)
                #logger.info("iqde %s --> pcu %f", iqde, pcu)
                return pcu
            res = so.minimize(
                self.culoss, startvals, method='SLSQP',  # trust-constr
                bounds=self.bounds,
                #            jac=gradient_respecting_bounds(self.bounds, self.culoss),
                constraints=[
                    {'type': 'eq',
                     'fun': lambda iqd: self.torque_iqd(*iqd) - torque}])
            #options={'disp': disp, 'maxiter': maxiter})
            if res['success']:
                return res.x
        logger.warning("%s: torque=%f %f, io=%s",
                       res['message'], torque, self.torque_iqd(*startvals),
                       startvals)
        raise ValueError(res['message'])

    def mtpa(self, i1max):
        """return iq, id, iex currents and maximum torque per current """
        T0 = self.torque_iqd(np.sqrt(2)*i1max, 0, self.bounds[-1][1])
        def i1tq(tq):
            return abs(i1max) - np.linalg.norm(self.iqd_torque(tq)[:2])/np.sqrt(2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tq = so.fsolve(i1tq, T0)[0]
            iq, id, iex = self.iqd_torque(tq)
        return iq, id, iex, tq

    def mtpa_tmech(self, i1max, n):
        """return iq, id, iex currents and maximum torque per current """
        T0 = self.torque_iqd(np.sqrt(2)*i1max, 0, self.bounds[-1][1])
        def i1tq(tq):
            return i1max - np.linalg.norm(self.iqd_tmech(tq, n)[:2])/np.sqrt(2)
        tq = so.fsolve(i1tq, T0)[0]
        iq, id, iex = self.iqd_tmech(tq, n)
        return iq, id, iex, tq

    def iqd_tmech_umax(self, torque, w1, u1max, log=0, **kwargs):
        """return currents and shaft torque at stator frequency and
         with minimal losses at max voltage"""
        iqde = self.iqd_tmech(torque, w1/2/np.pi/self.p)
        if np.linalg.norm(
                self.uqd(w1, *iqde)) <= u1max*np.sqrt(2):
            if log:
                log(iqde)
            return (*iqde, torque)
        #beta, i1 = betai1(iqde[0], iqde[1])
        #iex = iqde[2]

        #beta = 0 if torque>0 else np.pi
        io = iqde[0], 0, iqde[2] #*iqd(beta, i1), iex

        #    logger.debug("--- torque %g io %s", torque, io)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            n = w1/2/np.pi/self.p
            def sqrtculoss(iqde):
                pcu = self.culoss(iqde)
                #logger.info("iqde %s pcu %g", iqde, pcu)
                return pcu

            res = so.minimize(
                self.culoss, io, method='SLSQP',  # trust-constr
                bounds=self.bounds,
                constraints=[
                    {'type': 'eq',
                     'fun': lambda iqd: torque - self.tmech_iqd(*iqd, n)},
                    {'type': 'eq',
                     'fun': lambda iqd: u1max*np.sqrt(2)
                       - np.linalg.norm(self.uqd(w1, *iqd))}])
            #if res['success']:
            if log:
                log(res.x)
            return *res.x, self.tmech_iqd(*res.x, n)
        #logger.warning("%s: w1=%f torque=%f, u1max=%f, io=%s",
        #               res['message'], w1, torque, u1max, io)
        #raise ValueError(res['message'])
        #return [float('nan')]*4

    def iqd_torque_umax(self, torque, w1, u1max,
                        disp=False, maxiter=500, log=0, **kwargs):
        """return currents for torque with minimal losses"""
        iqde = self.iqd_torque(torque, disp, maxiter)
        if np.linalg.norm(
                self.uqd(w1, *iqde)) <= u1max*np.sqrt(2):
                if log:
                    log(iqde)
                return (*iqde, torque)
        io = iqde[0], 0, iqde[2]
        #    logger.debug("--- torque %g io %s", torque, io)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            def sqrtculoss(iqde):
                pcu = self.culoss(iqde)
                #logger.info("iqde %s pcu %g", iqde, pcu)
                return pcu

            res = so.minimize(
                self.culoss, io, method='SLSQP',  # trust-constr
                bounds=self.bounds,
                #options={'disp': disp, 'maxiter': maxiter},
                #            jac=gradient_respecting_bounds(self.bounds, self.culoss),
                constraints=[
                    {'type': 'eq',
                     'fun': lambda iqd: self.torque_iqd(*iqd) - torque},
                    {'type': 'eq',
                     'fun': lambda iqd: u1max*np.sqrt(2) - np.linalg.norm(
                         self.uqd(w1, *iqd))}])
            #if res['success']:
            if log:
                log(res.x)
            return *res.x, self.torque_iqd(*res.x)
            #logger.warning("%s: w1=%f torque=%f, u1max=%f, io=%s",
            #               res['message'], w1, torque, u1max, io)
            #raise ValueError(res['message'])

    def w1_imax_umax(self, i1max, u1max):
        """return frequency w1 and shaft torque at voltage u1max and current i1max

        Keyword arguments:
            u1max -- the maximum voltage (Vrms)
            i1max -- the maximum current (Arms)"""
        iq, id, iex, T = self.mtpa(i1max)
        n0 = np.sqrt(2)*u1max/np.linalg.norm(
            self.psi(iq, id, iex))/2/np.pi/self.p
        return self.w1_umax(u1max, iq, id, iex), T

    def w1_umax(self, u, iq, id, iex):
        """return frequency w1 at given voltage u and id, iq current

        Keyword arguments:
            u -- the maximum voltage (RMS)
            iq, id -- the d-q currents"""
        w10 = np.sqrt(2)*u/np.linalg.norm(self.psi(iq, id, iex))
        return so.fsolve(
            lambda w1: np.linalg.norm(self.uqd(w1, iq, id, iex))-u*np.sqrt(2),
            w10)[0]

    def characteristics(self, T, n, u1max, nsamples=50,
                        with_tmech=True, with_torque_corr=False, **kwargs):
        """calculate torque speed characteristics.
        return dict with list values of
        n, T, u1, i1, beta, cosphi, pmech, n_type

        Keyword arguments:
        T -- (float) the maximum torque in Nm
        n -- (float) the maximum speed in 1/s
        u1max -- (float) the maximum voltage in V rms
        nsamples -- (optional) number of speed samples
        with_tmech -- (optional) use friction and windage losses
        with_torque_corr -- (optional) T is corrected if out of range
        Optional arguments:
        i1max: max. phase current (RMS)
        """
        if kwargs.get('i1max', 0):
            w1type, T = self.w1_imax_umax(kwargs['i1max'], u1max)
        try:
            iq, id, iex = self.iqd_torque(T)
        except ValueError:
            tmax = self.torquemax(
                    self.i1range[1], self.exc_max)
            tmin = 0
            if self.betarange[0] < -np.pi/2:
                tmin = -self.torquemin(
                        self.i1range[1], self.exc_max)
            if with_torque_corr:
                Torig = T
                if T > 0:
                    T = np.floor(0.94*tmax)
                else:
                    T = np.ceil(0.94*tmin)
                logger.warning("corrected torque %f -> %f Nm",
                               Torig, T)
                iq, id, iex = self.iqd_torque(T)
            else:
                raise ValueError(
                        f"torque {T} Nm out of range ({tmin:.1f}, {tmax:.1f} Nm)")

        if with_tmech:
            i1max = betai1(iq, id)[1]
            if T < 0:
                i1max = -i1max
            w1type, Tf = self.w1_imax_umax(i1max, u1max)

        else:
            Tf = T
            w1type = self.w1_umax(u1max, iq, id, iex)
        logger.debug("w1type %f", w1type)
        wmType = w1type/self.p
        pmax = Tf*wmType

        def tload(wm):
            if abs(wm*Tf) < abs(pmax):
                return Tf
            return pmax/wm

        wmtab = []
        dw = 0
        wmMax = 5*wmType
        if n > 0:
            wmMax = min(wmMax, 2*np.pi*n)

        if wmType > wmMax:
            wmrange = sorted([0, wmMax])
            wmtab = np.linspace(0, wmMax, nsamples).tolist()
        else:
            wmrange = sorted([0, wmType, wmMax])
            nx = int(round(nsamples*wmType/wmMax))
            dw = (wmMax-wmType)/(nsamples-nx)
            wmtab = (np.linspace(0, wmType, nx).tolist() +
                     np.linspace(wmType+dw, wmMax, nsamples-nx).tolist())

        logger.info("Speed range T %f %s", Tf, wmrange)
        wmtab[0] = 0

        r = dict(u1=[], i1=[], id=[], iq=[], iex=[], T=[], cosphi=[], n=[],
                 beta=[], plfe1=[], plfe2=[], plcu1=[], plcu2=[])
        # add type speed to result dict
        r['n_type'] = wmType/2/np.pi
        for wm, tq in zip(wmtab, [tload(wx) for wx in wmtab]):
            w1 = wm*self.p
            if with_tmech:
                iq, id, iex, tqx = self.iqd_tmech_umax(
                        tq, w1, u1max)
            else:
                iq, id, iex, tqx = self.iqd_torque_umax(
                        tq, w1, u1max)
                tqx -= self.tfric
            uq, ud = self.uqd(w1, iq, id, iex)
            u1 = np.linalg.norm((uq, ud))/np.sqrt(2)
            f1 = w1/2/np.pi
            r['id'].append(id)
            r['iq'].append(iq)
            r['iex'].append(iex)
            r['u1'].append(u1)
            beta = np.arctan2(id, iq)
            if beta > 0:
                beta -= 2*np.pi
            i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
            r['i1'].append(i1)
            r['beta'].append(beta/180*np.pi)
            gamma = np.arctan2(ud, uq)
            r['cosphi'].append(np.cos(gamma - beta))
            r['plfe1'].append(self.iqd_plfe1(iq, id, iex, f1))
            r['plfe2'].append(self.iqd_plfe2(iq, id, iex, f1))
            r['plcu1'].append(self.m*i1**2*self.rstat(w1))
            r['plcu2'].append(iex**2*self.rrot(0))
            r['T'].append(tqx)
            r['n'].append(wm/2/np.pi)

        r['plfe'] = (np.array(r['plfe1']) + np.array(r['plfe2'])).tolist()
        r['plcu'] = (np.array(r['plcu1']) + np.array(r['plcu2'])).tolist()
        r['plfw'] = [self.pfric(n) for n in r['n']]
        r['pmech'] = [2*np.pi*n*tq
                      for n, tq in zip(r['n'], r['T'])]
        pmech = np.array(r['pmech'])
        pltotal = (np.array(r['plfe1']) + np.array(r['plfw']) +
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

class SynchronousMachinePsidq(SynchronousMachine):

    def __init__(self, eecpars, lfe=1, wdg=1, **kwargs):
        super(self.__class__, self).__init__(
            eecpars, **kwargs)
        self.iqrange = (eecpars['psidq'][0]['iq'][0],
                        eecpars['psidq'][0]['iq'][-1])
        self.idrange = (eecpars['psidq'][0]['id'][0],
                        eecpars['psidq'][0]['id'][-1])

        self.betarange = (-np.pi if min(self.iqrange) < 0 else -np.pi/2,
                          0 if max(self.iqrange) > 0 else -np.pi/2)
        self.i1range = (0, betai1(np.max(self.iqrange), 0)[1])
        islinear = True
        iexc = [l['ex_current'] for l in eecpars['psidq']]
        self.exc_max = iexc[-1]
        id = [i/wdg for i in eecpars['psidq'][-1]['id']]
        idx = np.linspace(id[0], id[-1], 12)
        iq = [i/wdg for i in eecpars['psidq'][-1]['iq']]
        iqx = np.linspace(iq[0], iq[-1], 20)

        if _islinear(iexc):
            exc = iexc
            psid = wdg*lfe*np.array([
                _splinterp(iq, id, iqx, idx, l['psid'])
                for l in eecpars['psidq']])
            psiq = wdg*lfe*np.array([
                _splinterp(iq, id, iqx, idx, l['psiq'])
                for l in eecpars['psidq']])
        else:
            islinear = False
            nsamples = 10
            iexcl = np.linspace(iexc[0], iexc[-1], nsamples)
            exc = iexcl
            psid = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                [_splinterp(iq, id, iqx, idx, l['psid'])
                 for l in eecpars['psidq']]))
            psiq = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                [_splinterp(iq, id, iqx, idx, l['psiq'])
                 for l in eecpars['psidq']]))

        self.psidf = ip.RegularGridInterpolator(
            (exc, iqx, idx), psid,
            method='cubic', bounds_error=False, fill_value=None)
        self.psiqf = ip.RegularGridInterpolator(
            (exc, iqx, idx), psiq,
            method='cubic', bounds_error=False, fill_value=None)
        self.bounds = [(min(iq), max(iq)),
                       (min(id), 0),
                       (iexc[0], iexc[-1])]
        # iron losses
        idname = 'psidq'
        keys = [k for k in self.plexp.keys() if k in eecpars[idname][0]['losses']]
        try:
            # check if bertotti
            if 'styoke_excess' in eecpars[idname][0]['losses'] and \
               np.any(np.array(eecpars[idname][0]['losses']['styoke_excess'])):
                self.bertotti = True
                keys += ['styoke_excess', 'stteeth_excess','rotor_excess']
            if islinear:
                pfe = {k: np.array([l['losses'][k]
                                    for l in eecpars[idname]])
                       for k in keys}
            else:
                pfe = {k: _linsampl(iexc, iexcl,
                                    np.array([l['losses'][k]
                                             for l in eecpars[idname]]))
                       for k in keys}
                self._losses = {k: ip.RegularGridInterpolator(
                    (exc, iq, id), lfe*pfe[k],
                    method='cubic', bounds_error=False, fill_value=None)
                            for k in keys}
                self._set_losspar(eecpars[idname][0]['losses']['speed'],
                                  eecpars[idname][0]['losses']['ef'],
                                  eecpars[idname][0]['losses']['hf'])
        except KeyError:
            logger.warning("loss map missing")
            self._losses = {k: lambda x: 0 for k in (
                'styoke_hyst', 'stteeth_hyst',
                'styoke_eddy', 'stteeth_eddy',
                'rotor_hyst', 'rotor_eddy')}


    def psi(self, iq, id, iex):
        """return psid, psiq of currents iq, id"""
        try:
            return self.psidf((iex, iq, id)), self.psiqf((iex, iq, id))
        except ValueError as ex:
            logger.error(iex, iq, id)
            raise ex

    def plfe1(self, iq, id, iex, f1):
        losskeys = [k for k in self._losses if k in [
            'styoke_eddy', 'styoke_hyst',
            'stteeth_eddy', 'stteeth_hyst']]
        if self.bertotti:
            losskeys += ['styoke_excess', 'stteeth_excess']
        return np.sum([
            self._losses[k]((iex, iq, id))*(f1/self.fo)**self.plexp[k][0]
            for k in losskeys], axis=0)

    def plfe2(self, iq, id, iex, f1):
        losskeys = ['rotor_hyst', 'rotor_eddy']
        if self.bertotti:
            losskeys += ['rotor_excess']
        return np.sum([
            self._losses[k]((iex, iq, id))*(f1/self.fo)**self.plexp[k][0]
            for k in losskeys], axis=0)

    def iqd_plfe1(self, iq, id, iex, f1):
        return self.plfe1(iq, id, iex, f1)

    def iqd_plfe2(self, iq, id, iex, f1):
        return self.plfe2(iq, id, iex, f1)

class SynchronousMachineLdq(SynchronousMachine):
    def __init__(self, eecpars, lfe=1, wdg=1, **kwargs):
        super(self.__class__, self).__init__(eecpars, **kwargs)
        self.betarange = (eecpars['ldq'][0]['beta'][0]/180*np.pi,
                          eecpars['ldq'][0]['beta'][-1]/180*np.pi)
        self.i1range = (0, eecpars['ldq'][0]['i1'][-1])
        islinear = True
        iexc = [l['ex_current'] for l in eecpars['ldq']]
        self.exc_max = iexc[-1]
        i1 = [i/wdg for i in eecpars['ldq'][-1]['i1']]
        i1x = np.linspace(i1[0], i1[-1], 12)
        beta = [b*np.pi/180 for b in eecpars['ldq'][-1]['beta']]
        betax = np.linspace(beta[0], beta[-1], 20)

        if _islinear(iexc):
            logger.info("Linear sampled ex current: %s",
                        iexc)
            exc = iexc
            psid = wdg*lfe*np.array([
                _splinterp(beta, i1, betax, i1x, l['psid'])
                    for l in eecpars['ldq']])
            psiq = wdg*lfe*np.array([
                _splinterp(beta, i1, betax, i1x, l['psiq'])
                    for l in eecpars['ldq']])
        else:
            islinear = False
            logger.info("Non Linear sampled ex current %s",
                        iexc)
            nsamples = 10
            iexcl = np.linspace(iexc[0], iexc[-1], nsamples)
            exc = iexcl
            psid = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                [_splinterp(beta, i1, betax, i1x, l['psid'])
                     for l in eecpars['ldq']]))
            psiq = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                [_splinterp(beta, i1, betax, i1x, l['psiq'])
                     for l in eecpars['ldq']]))

        # extrapolate outside range
        self.psidf = ip.RegularGridInterpolator(
            (exc, betax, i1x), np.sqrt(2)*psid,
            method='cubic',
            bounds_error=False, fill_value=None)
        self.psiqf = ip.RegularGridInterpolator(
            (exc, betax, i1x), np.sqrt(2)*psiq,
            method='cubic'
            , bounds_error=False, fill_value=None)
        i1max = np.sqrt(2)*(max(i1))
        self.bounds = [(np.cos(min(beta))*i1max, i1max),
                       (-i1max, 0),
                       (iexc[0], iexc[-1])]

        # iron losses
        idname = 'ldq'
        keys = [k for k in self.plexp.keys() if k in eecpars[idname][0]['losses']]
        try:
            # check bertotti losses
            if 'styoke_excess' in eecpars[idname][0]['losses'] and \
               np.any(np.array(eecpars[idname][0]['losses']['styoke_excess'])):
                self.bertotti = True
                keys += ['styoke_excess', 'stteeth_excess','rotor_excess']

            if islinear:
                pfe = {k: np.array([l['losses'][k]
                                    for l in eecpars[idname]])
                       for k in keys}
            else:
                pfe = {k: _linsampl(iexc, iexcl,
                                    np.array([l['losses'][k]
                                             for l in eecpars[idname]]))
                       for k in keys}

        # fill value with nan outside range
            self._losses = {k: ip.RegularGridInterpolator(
                (exc, beta, i1), lfe*pfe[k],
                method='cubic', bounds_error=False, fill_value=None)
                            for k in pfe.keys()}
            self._set_losspar(eecpars[idname][0]['losses']['speed'],
                              eecpars[idname][0]['losses']['ef'],
                              eecpars[idname][0]['losses']['hf'])
        except KeyError:
            logger.warning("loss map missing")
            self._losses = {k: lambda x: 0 for k in (
                'styoke_hyst', 'stteeth_hyst',
                'styoke_eddy', 'stteeth_eddy',
                'rotor_hyst', 'rotor_eddy')}
            pass

    def psi(self, iq, id, iex):
        """return psid, psiq of currents iq, id"""
        beta = np.arctan2(id, iq)
        if beta > 0:
            beta -= 2*np.pi
        i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
        try:
            return self.psidf((iex, beta, i1)), self.psiqf((iex, beta, i1))
        except ValueError as ex:
            logger.error("iex %s iq %f id %f beta %f i1 %f",
                         iex, iq, id, beta, i1)
            raise ex

    def plfe1(self, beta, i1, iex, f1):
        losskeys = ['styoke_eddy', 'styoke_hyst',
                     'stteeth_eddy', 'stteeth_hyst']
        if self.bertotti:
            losskeys += ['styoke_excess', 'stteeth_excess']
        return np.sum([
            self._losses[k]((iex, beta, i1))*(f1/self.fo)**self.plexp[k][0]
            for k in losskeys if k in self._losses], axis=0)

    def plfe2(self, beta, i1, iex, f1):
        losskeys =  ['rotor_hyst', 'rotor_eddy']
        if self.bertotti:
            losskeys += ['rotor_excess']
        return np.sum([
            self._losses[k]((iex, beta, i1))*(f1/self.fo)**self.plexp[k][0]
            for k in losskeys], axis=0)

    def iqd_plfe1(self, iq, id, iex, f1):
        beta = np.arctan2(id, iq)
        if np.isscalar(beta):
            if beta > 0:
                beta -= 2*np.pi
        else:
            beta[beta > 0] -= 2*np.pi
        i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
        return self.plfe1(beta, i1, iex, f1)

    def iqd_plfe2(self, iq, id, iex, f1):
        beta = np.arctan2(id, iq)
        if np.isscalar(beta):
            if beta > 0:
                beta -= 2*np.pi
        else:
            beta[beta > 0] -= 2*np.pi
        i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
        return self.plfe2(beta, i1, iex, f1)


if __name__ == '__main__':
    import sys
    import json
    import femagtools.plot
    import matplotlib.pyplot as plt
    with open(sys.argv[1]) as fp:
        eecpar = json.load(fp)
    if 'ldq' in eecpar:
        m = SynchronousMachineLdq(eecpar)
    else:
        m = SynchronousMachinePsidq(eecpar)
        T = 240
        u1max = 163
        nmax = 1000
        r = m.characteristics(T, 0, u1max)
        femagtools.plot.characteristics(r, 'SynchronousMachine')
        plt.show()
