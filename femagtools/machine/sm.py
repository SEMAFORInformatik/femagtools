"""
  femagtools.sm
  ~~~~~~~~~~~~~
  wound-rotor synchronous machine (EESM) electrical circuit model

  Copyright 2022: Semafor Informatik & Energie AG, Switzerland
"""
import logging
import warnings
import numpy as np
import scipy.optimize as so
import scipy.interpolate as ip
from .utils import skin_resistance, wdg_resistance
from .. import parstudy, windings
import femagtools.bch

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


def parident(workdir, engine, machine,
             magnetizingCurves, condMat,
             **kwargs):
        """return dict of parameters of equivalent circuit for
        electrically excited synchronous machines

        arguments:
        machine -- dict() with machine parameters
        magnetizingCurves -- list of dict() with BH curves
        condMat -- list of dict() with conductor material properties

        optional arguments:
        num_cur_steps: number of current steps (default 5)
        num_beta_steps: number of current steps (default 13)
        num_exc_steps: number of excitation current (default 7)
        speed: rotor speed in 1/s (default 160/p)
        i1_max: maximum current in A rms (default approx 3*i1nom)
        """
        da1 = machine['outer_diam']
        Q1 = machine['stator']['num_slots']
        if 'statorRotor3' in machine['stator']:
            hs = machine['stator']['statorRotor3']['slot_height']
        elif 'stator1' in machine['stator']:
            hs = machine['stator']['stator1']['slot_rf1'] - machine['stator']['stator1']['tip_rh1']
        elif 'stator4' in machine['stator']:
            hs = machine['stator']['stator4']['slot_height']
        N = machine['windings']['num_wires']
        Jmax = 15  # max current density in A/mm2

        i1_max = round(0.28*np.pi*hs*(da1+hs)/Q1/N*Jmax*1e5)*10 * \
            machine['windings'].get('num_par_wdgs', 1)

        ifnom = machine['rotor']['ifnom']
        exc_logspace = True
        if exc_logspace:
            excur = np.logspace(np.log(ifnom/10), np.log(1.5*ifnom),
                                kwargs.get("num_exc_steps", 8),
                                base=np.exp(1)).tolist()
        else:
            excur = np.linspace(ifnom/10, 1.5*ifnom,
                                kwargs.get("num_exc_steps", 8))

        parvardef = {
            "decision_vars": [
                {"values": excur, "name": "load_ex_cur"}
            ]
        }

        parvar = parstudy.List(
            workdir,  condMat=condMat,
            magnetizingCurves=magnetizingCurves)

        simulation = dict(
                calculationMode=kwargs.get('calculationMode',
                                           'ld_lq_fast'),
                wind_temp=20.0,
                i1_max=kwargs.get('i1_max', i1_max),
                maxid=0,
                minid=-i1_max,
                maxiq=i1_max,
                miniq=-i1_max,
                delta_id=i1_max/kwargs.get('num_cur_steps', 5),
                delta_iq=i1_max/kwargs.get('num_cur_steps', 5),
                beta_min=-180.0,
                beta_max=0.0,
                calc_noload=0,
                num_move_steps=kwargs.get('num_move_steps', 31),
                load_ex_cur=0.5,
                num_cur_steps=kwargs.get('num_cur_steps', 5),
                num_beta_steps=kwargs.get('num_beta_steps', 13),
                num_par_wdgs=machine['windings'].get('num_par_wdgs', 1),
                period_frac=6,
                speed=50.0)

        ###self.cleanup()  # remove previously created files in workdir
        results = parvar(parvardef, machine, simulation, engine)
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
        b = results['f'][-1]

        losskeys = ('speed', 'ef', 'hf',
                    'styoke_hyst', 'stteeth_hyst', 'styoke_eddy',
                    'stteeth_eddy', 'rotor_hyst', 'rotor_eddy')

        # winding resistance
        try:
            r1 = machine['windings']['resistance']
        except KeyError:
            yd = machine['windings'].get('coil_span', Q1/machine['poles'])
            wdg = windings.Winding(
            {'Q': machine['stator']['num_slots'],
             'm': machine['windings']['num_phases'],
             'p': machine['poles']//2,
             'l': machine['windings']['num_layers'],
             'yd': yd})

            lfe = machine['lfe']
            g = machine['windings'].get('num_par_wdgs', 1)
            if 'dia_wire' in machine['windings']:
                aw = np.pi*machine['windings'].get('dia_wire', 1e-3)**2/4
            else:  # wire diameter from slot area
                aw = 0.75 * \
                        machine['windings'].get('cufilfact', 0.45)*np.pi*da1*hs/Q1/2/N
            r1 = wdg_resistance(wdg, N, g, aw, da1, hs, lfe)

        if simulation['calculationMode'] == 'ld_lq_fast':
            return dict(m=3, p=b['machine']['p'],
                        r1=r1,
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
                                        for k in losskeys})
                             for b in results['f']])
        return dict(m=3, p=b['machine']['p'],
                    r1=r1,
                    r2=machine['rotor'].get('resistance', 1),
                    rotor_mass=rotor_mass, kfric_b=1,
                    psidq=[dict(
                            ex_current=b['machine']['ex_current'],
                            iq=b['psidq']['iq'],
                            id=b['psidq']['id'],
                            psid=b['psidq']['psid'],
                            psiq=b['psidq']['psiq'],
                            torque=b['psidq']['torque'],
                            losses={k: b['psidq']['losses'][k]
                                    for k in losskeys})
                           for b in results['f']])

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
    def __init__(self, eecpars):
        for k in eecdefaults.keys():
            setattr(self, k, eecdefaults[k])

        for k in eecpars:
            if k not in ('ldq', 'psidq'):
                setattr(self, k, eecpars[k])

        self.fo = 50
        self.plexp = {'styoke_hyst': 1.0,
                      'stteeth_hyst': 1.0,
                      'styoke_eddy': 2.0,
                      'stteeth_eddy': 2.0,
                      'rotor_hyst': 1.0,
                      'rotor_eddy': 2.0}

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
        return self.plcu1(iqde, w1) + self.plcu2(iqde, 0)

    def iqd_plcu1(self, iq, id, f1):
        return self.plcu1((iq, id), 2*np.pi*f1)

    def iqd_plcu2(self, iq, id, iex):
        return self.plcu2((iq, id, iex))

    def iqd_plfe2(self, iq, id, f1):
        return np.zeros(np.asarray(iq).shape)

    def iqd_plmag(self, iq, id, f1):
        return np.zeros(np.asarray(iq).shape)

    def iqd_torque(self, torque, disp=False, maxiter=500):
        """return currents for torque with minimal losses"""
        if torque > 0:
            startvals = self.bounds[0][1]/2, 0, sum(self.bounds[-1])/2
        else:
            startvals = -self.bounds[0][1]/2, 0, sum(self.bounds[-1])/2

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            def sqrtculoss(iqde):
                pcu = self.culoss(iqde)
                #logger.info("iqde %s --> pcu %f", iqde, pcu)
                return pcu
            res = so.minimize(
                sqrtculoss, startvals, method='SLSQP',  # trust-constr
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
        #logger.info(">> torque %g w1 %g u1 %g io %s", torque, w1, u1max, io)
        if io == 0:
            iqde = self.iqd_torque(torque, disp, maxiter)
            if np.linalg.norm(
                    self.uqd(w1, *iqde)) <= u1max*np.sqrt(2):
                return (*iqde, torque)
            io = iqde[0], 0, iqde[2]
            logger.debug("--- torque %g io %s", torque, io)
        #logger.info(">>      io %s", io)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            def sqrtculoss(iqde):
                pcu = self.culoss(iqde)
                #logger.info("iqde %s pcu %g", iqde, pcu)
                return pcu

            res = so.minimize(
                sqrtculoss, io, method='SLSQP',  # trust-constr
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

        logger.info("Speed range %s", wmrange)
        wmtab[0] = 0

        r = dict(u1=[], i1=[], id=[], iq=[], iex=[], T=[], cosphi=[], n=[],
                 beta=[], plfe1=[], plcu1=[], plcu2=[])
        T = [tload(wx) for wx in wmtab]
        tfric = self.kfric_b*self.rotor_mass*30e-3/np.pi
        w1tab = []
        for wm, tq in zip(wmtab, T):
            #            try:
            w1 = wm*self.p
            tqx = tq
            #            if w1 <= w1type:
            #                iq, id, iex = self.iqd_torque(tq)
            #            else:
            iq, id, iex, tqx = self.iqd_torque_umax(
                tq, w1, u1max)
                #                        (0.9*iq, 0.9*id,
                #                         min(self.bounds[-1][0], 0.9*iex)))[:-1]
            #logger.info("w1 %g tq %g: iq %g iex %g tqx %g",
            #            w1, tq, iq, iex, tqx)
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
            r['plcu1'].append(self.m*i1**2*self.rstat(w1))
            r['plcu2'].append(iex**2*self.rrot(0))
            r['T'].append(tq-tfric)
            r['n'].append(wm/2/np.pi)
            # except ValueError as ex:
            #    logger.warning("ex %s wm %f T %f", ex, wm, tq)
            #    break

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


class SynchronousMachinePsidq(SynchronousMachine):

    def __init__(self, eecpars, lfe=1, wdg=1):
        super(self.__class__, self).__init__(
                eecpars)
        self.iqrange = (eecpars['psidq'][0]['iq'][0],
                        eecpars['psidq'][0]['iq'][-1])
        self.idrange = (eecpars['psidq'][0]['id'][0],
                        eecpars['psidq'][0]['id'][-1])
        islinear = True
        iexc = [l['ex_current'] for l in eecpars['psidq']]
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
        try:
            idname = 'psidq'
            keys = self.plexp.keys()
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
                    (exc, iq, id), lfe*np.array(pfe[k]),
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
        return np.sum([
            self._losses[k]((iex, iq, id))*(f1/self.fo)**self.plexp[k][0]
            for k in ('styoke_eddy', 'styoke_hyst',
                      'stteeth_eddy', 'stteeth_hyst')], axis=0)

    def iqd_plfe1(self, iq, id, iex, f1):
        return self.plfe1(iq, id, iex, f1)


class SynchronousMachineLdq(SynchronousMachine):
    def __init__(self, eecpars, lfe=1, wdg=1):
        super(self.__class__, self).__init__(eecpars)
        self.betarange = (eecpars['ldq'][0]['beta'][0]/180*np.pi,
                          eecpars['ldq'][0]['beta'][-1]/180*np.pi)
        self.i1range = (0, eecpars['ldq'][0]['i1'][-1])
        islinear = True
        iexc = [l['ex_current'] for l in eecpars['ldq']]
        i1 = [i/wdg for i in eecpars['ldq'][-1]['i1']]
        i1x = np.linspace(i1[0], i1[-1], 12)
        beta = [b*np.pi/180 for b in eecpars['ldq'][-1]['beta']]
        betax = np.linspace(beta[0], beta[-1], 20)

        if _islinear(iexc):
            exc = iexc
            psid = wdg*lfe*np.array([
                    _splinterp(beta, i1, betax, i1x, l['psid'])
                    for l in eecpars['ldq']])
            psiq = wdg*lfe*np.array([
                    _splinterp(beta, i1, betax, i1x, l['psiq'])
                    for l in eecpars['ldq']])
        else:
            islinear = False
            nsamples = 10
            iexcl = np.linspace(iexc[0], iexc[-1], nsamples)
            exc = iexcl
            psid = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                    [_splinterp(beta, i1, betax, i1x, l['psid'])
                     for l in eecpars['ldq']]))
            psiq = wdg*lfe*_linsampl(iexc, iexcl, np.array(
                    [_splinterp(beta, i1, betax, i1x, l['psiq'])
                     for l in eecpars['ldq']]))

        self.psidf = ip.RegularGridInterpolator(
                (exc, betax, i1x), np.sqrt(2)*psid,
                method='cubic', bounds_error=False, fill_value=None)
        self.psiqf = ip.RegularGridInterpolator(
                (exc, betax, i1x), np.sqrt(2)*psiq,
                method='cubic', bounds_error=False, fill_value=None)
        i1max = np.sqrt(2)*(max(i1))
        self.bounds = [(np.cos(min(beta))*i1max, i1max),
                       (-i1max, 0),
                       (iexc[0], iexc[-1])]
        keys = self.plexp.keys()
        try:
            idname = 'ldq'
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
                    (exc, beta, i1), lfe*np.array(pfe[k]),
                    method='cubic', bounds_error=False, fill_value=None)
                            for k in keys}
            self._set_losspar(eecpars[idname][0]['losses']['speed'],
                              eecpars[idname][0]['losses']['ef'],
                              eecpars[idname][0]['losses']['hf'])
        except ValueError:
            logger.warning("loss map missing")
            self._losses = {k: lambda x: 0 for k in (
                'styoke_hyst', 'stteeth_hyst',
                'styoke_eddy', 'stteeth_eddy',
                'rotor_hyst', 'rotor_eddy')}

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

    def plfe1(self, beta, i1, iex, f1):
        return np.sum([
            self._losses[k]((iex, beta, i1))*(f1/self.fo)**self.plexp[k][0]
            for k in ('styoke_eddy', 'styoke_hyst',
                      'stteeth_eddy', 'stteeth_hyst')], axis=0)
    def iqd_plfe1(self, iq, id, iex, f1):
        beta = np.arctan2(id, iq)
        if np.isscalar(beta):
            if beta > 0:
                beta -= 2*np.pi
        else:
            beta[beta > 0] -= 2*np.pi
        i1 = np.linalg.norm((id, iq), axis=0)/np.sqrt(2.0)
        return self.plfe1(beta, i1, iex, f1)


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
