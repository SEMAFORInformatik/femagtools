import pathlib
import logging
import json
import numpy as np
import scipy.optimize as so
from scipy.interpolate import make_interp_spline
import scipy.integrate as ig
import femagtools.parstudy

logger = logging.getLogger('shortcircuit')

def get_ik(nc, i1, Hk, eps=5e-4):
    """return approx. current at knee point in demag characteristics
    based on 2 operating points (noload, load)"""
    def find_xpeak(x, y):
        """ return first peak from the right
        """
        peaks = (np.diff(np.sign(np.diff(y))) < 0).nonzero()[0] + 1
        ip = np.where(y[peaks] > eps)[0][-1]
        return x[peaks][ip], y[peaks][ip]

    elements = nc.magnet_elements()
    icur = -1
    ibeta = -1
    ncdemag = np.array([nc.demagnetization(e, icur, ibeta)[1]
                        for e in elements])
    nbins = 50
    imax = np.argmax(np.max(ncdemag, axis=0))
    n, b = np.histogram(ncdemag[:, imax], bins=nbins, density=True)
    x, y = find_xpeak(b,n)
    icur = 0
    ibeta = 0
    ncdemag = np.array([nc.demagnetization(e, icur, ibeta)[1]
                        for e in elements])
    imax = np.argmax(np.max(ncdemag, axis=0))
    n, b = np.histogram(ncdemag[:, imax], bins=nbins, density=True)

    y0 = b[np.argmax(n)]/-Hk
    y1 = x/-Hk
    x1 = np.sqrt(2)*i1
    m = (y1-y0)/x1
    return (1 - y0)/m

def _parstudy_list(femag, result_func):
    workdir = femag.workdir
    magnetMat = femag.magnets
    magnetizingCurves = femag.magnetizingCurves
    condMat = femag.condMat
    templatedirs = femag.templatedirs
    cmd = femag.cmd
    return femagtools.parstudy.List(
        workdir, condMat=condMat, magnets=magnetMat,
        magnetizingCurves=magnetizingCurves,
        templatedirs=templatedirs,
        cmd=cmd, result_func=result_func)

def get_shortCircuit_parameters(bch, nload):
    """extracts shortcircuit parameters from bch"""
    try:
        if nload < 0:
            nload = 0
        if nload > 2:
            nload = 2
        if nload > 0:
            dqld = bch.dqPar['ld']
            dqlq = bch.dqPar['lq']
            dqpsim = bch.dqPar['psim']
            if len(dqld) <= nload or len(dqlq) <= nload or len(dqpsim) <= nload:
                ld = dqld[-1]/bch.armatureLength
                lq = dqlq[-1]/bch.armatureLength
                psim = dqpsim[-1]/bch.armatureLength
            else:
                ld = dqld[nload-1]/bch.armatureLength
                lq = dqlq[nload-1]/bch.armatureLength
                psim = dqpsim[nload-1]/bch.armatureLength
        else:
            ld = bch.machine['ld']/bch.armatureLength
            lq = bch.machine['lq']/bch.armatureLength
            psim = bch.machine['psim']/bch.armatureLength

        return dict(
            r1=bch.machine['r1'],
            l_endwinding=bch.machine.get('ls1', 0),
            ld=ld,
            lq=lq,
            Hk=bch.magnet['demag_hx'],
            Tmag=bch.magnet['Tmag'],
            psim=psim,
            current_angles=bch.current_angles,
            num_pol_pair=bch.machine['p'],
            fc_radius=bch.machine['fc_radius'],
            lfe=bch.armatureLength/1e3,
            pocfilename=bch.machine['pocfile'],
            num_par_wdgs=bch.machine.get('num_par_wdgs', 1),
            calculationMode='shortcircuit')
    except (KeyError, AttributeError, IndexError):
        raise ValueError("missing pm/Rel-Sim results")


def find_peaks_and_valleys(t, y):
    """ return peak and valley of y with maximum amplitude
    """
    peaks = (np.diff(np.sign(np.diff(y))) < 0).nonzero()[0] + 1
    if len(peaks>0):
        ip = np.argmax(y[peaks])
        pv = {'ip': y[peaks][ip], 'tp': t[peaks][ip]}
    else:
        pv = {'ip': [], 'tp': []}
    valleys = (np.diff(np.sign(np.diff(y))) > 0).nonzero()[0] + 1
    if len(valleys>0):
        iv = np.argmin(y[valleys])
        pv.update({'iv': y[valleys][iv], 'tv': t[valleys][iv]})
    else:
        pv.update({'iv': [], 'tv': []})
    pv.update({'peaks': y[peaks], 'valleys': y[valleys],
               'tpeaks': t[peaks], 'tvalleys': t[valleys]})
    return pv


def shortcircuit(femag, machine, bch, simulation, engine=0):
    "run a 2 or 3ph short circuit simulation"
    scdata = {}
    sccalcmode = simulation.get('calculationMode', '')
    simulation.update(
        get_shortCircuit_parameters(bch,
                                    simulation.get('initial', 2)))
    try:
        from .machine.utils import KTH
        r1 = femag.model.winding['resistance']
        simulation['r1'] = r1*(1+KTH*(simulation['wind_temp']-20))
    except (KeyError, AttributeError):
        pass
    if 'speed' not in simulation:
        simulation['speed'] = bch.dqPar['speed']

    if simulation.get('sim_demagn', 0):
        nc = femag.read_nc()
        simulation['ik'] = get_ik(
            nc, bch.machine['i1'], simulation['Hk'])

    phirot = bch.flux['1'][0]['displ']
    #if (len(phirot)-1) % 3 == 0:
    #    phirot = phirot[::3]
    #elif (len(phirot)-1) % 2 == 0:
    #    phirot = phirot[::2]

    if simulation.get('sc_type', 3) == 3:
        logger.info("3phase short circuit simulation")
        builder = femagtools.fsl.Builder(femag.templatedirs)
        femag.model.calc_fe_loss = 0
        fslcmds = (builder.open_model(femag.model) +
                   builder.create_shortcircuit(simulation))
#                   ['save_model("cont")'])
        fslfile = 'shortcircuit.fsl'
        (pathlib.Path(femag.workdir)/fslfile).write_text(
            '\n'.join(fslcmds),
            encoding='latin1', errors='ignore')
        femag.run(fslfile)  # options?
        bchfile = femag.get_bch_file(femag.modelname)
        if bchfile:
            bchsc = femagtools.bch.Reader()
            logger.info("Read BCH %s", bchfile)
            bchsc.read(pathlib.Path(bchfile).read_text(
                encoding='latin1', errors='ignore'))
            bchsc.scData['demag'] = bchsc.demag
            if simulation.get('sim_demagn', 0):
                dd = {'displ': [d['displ']
                                for d in bchsc.demag if 'displ' in d],
                      'H_max': [d['H_max']
                                for d in bchsc.demag if 'H_max' in d],
                      'H_av': [d['H_av']
                               for d in bchsc.demag if 'H_av' in d]}
                istat = np.array([bchsc.scData[k]
                                  for k in ('ia', 'ib', 'ic')])
                k = np.argmax(np.abs(istat))
                i1max = simulation.get('i1max',
                                       np.max(np.abs(istat[:, k])))
                bchsc.scData['demag'] = demag(
                    femag, machine, simulation,
                    i1max, phirot)
                bchsc.scData['demag'].update(dd)
            scdata = bchsc.scData
            #for w in bch.flux:
            #    try:
            #        bch.flux[w] += bchsc.flux[w]
            #        bch.flux_fft[w] += bchsc.flux_fft[w]
            #    except (KeyError, IndexError):
            #        logging.debug(
            #            "No additional flux data in sc simulation")
            #        break

    if simulation.get('sc_type', 3) == 2:
        if (len(phirot)-1) % 3 == 0:
            phirot = phirot[::3]
        elif (len(phirot)-1) % 2 == 0:
            phirot = phirot[::2]
        simulation['phi'] = phirot
        if 'i1max' not in simulation:
            # just a wild guess
            simulation['i1max'] = 5*bch.machine['i1']
        logger.info("2phase short circuit simulation i1max = %.0f",
                    simulation['i1max'])
        if 'magn_temp' not in simulation:
            simulation['magn_temp'] = bch.magnet.get('Tmag', 20)
        scdata = shortcircuit_2phase(femag, machine, simulation, engine)

    else:
        logger.warning("Empty shortcircuit results for type %s",
                       simulation.get('sc_type', 'unknown'))
    # must reset sccalcmode
    if sccalcmode:
        simulation['calculationMode'] = sccalcmode
    else:
        del simulation['calculationMode']
    return scdata

def sc_result_func(task):
    basedir = pathlib.Path(task.directory)
    psitorq = np.loadtxt(basedir/'psi-torq-rot.dat')
    pos = np.unique(psitorq[:, 0])
    ncurs = psitorq[:, 0].shape[0]//pos.shape[0]
    ire = psitorq[:ncurs, 1:4]
    psire = np.reshape(psitorq[:, 4:7], (-1, ncurs, 3))
    torq = np.reshape(psitorq[:, 7], (-1, ncurs))
    pv = find_peaks_and_valleys(pos, psire[:, 0, 0])
    pmod = pv['peaks'].shape[0]
    return {'pos': pos.tolist(), 'ire': ire.tolist(),
            'psire': psire.tolist(), 'torq': torq.tolist(),
            'pmod': pmod}


def shortcircuit_2phase(femag, machine, simulation, engine=0):
    i1max = simulation['i1max']
    num_cur_steps = 4
    i1 = np.linspace(0, i1max, num_cur_steps)
    i1vec = np.concat((-i1[::-1], i1[1:]))
    simulation.update({
        'calculationMode': 'psi-torq-rot',
        'curvec': [],
        'magntemp': simulation['magn_temp']})

    if engine:
        parstudy = _parstudy_list(femag, sc_result_func)
        parvardef = {
            "decision_vars": [
                {"values": i1vec, "name": "curvec"}]
        }
        results = parstudy(parvardef, machine, simulation, engine)
        ire = np.array([r['ire'][0] for r in results['f']])
        pos = np.array(results['f'][0]['pos'])
        #phi = pos*np.pi/180
        torq = np.hstack([r['torq'] for r in results['f']])
        psire = np.hstack([r['psire'] for r in results['f']])
        simulation['curvec'] = i1vec.tolist()
    else:
        simulation['curvec'] = i1vec.tolist()
        results = femag(machine, simulation, fslfile='2ph_shortcircuit.fsl')
        class Task:
            def __init__(self, workdir):
                self.directory = workdir
        results = sc_result_func(Task(femag.workdir))
        ire = np.array(results['ire'])
        pos = np.array(results['pos'])
        torq = np.array(results['torq'])
        psire = np.array(results['psire'])

    pmod = 2
    #with open('results.json', 'w') as fp:
    #    json.dump({'ire': ire.tolist(), 'pos': pos.tolist(),
    #               'torq': torq.tolist(), 'psire': psire.tolist()}, fp)
    logger.info("move steps %d currents %s pmod %s",
                len(pos), ire[:,0], pmod)

    Ai = [femagtools.utils.fft(pos, psire[:, k, 0], pmod=pmod)['a']
          for k in range(np.shape(psire)[1])]
    A = make_interp_spline(ire[:,0], Ai)
    A0i = [femagtools.utils.fft(pos, psire[:, k, 0], pmod=pmod)['a0']
           for k in range(np.shape(psire)[1])]
    A0 = make_interp_spline(ire[:,0], A0i)
    Bi = [femagtools.utils.fft(pos, psire[:, k, 1], pmod=pmod)['a']
          for k in range(np.shape(psire)[1]-1, -1, -1)]
    B = make_interp_spline(ire[::-1,1], Bi)
    B0i = [femagtools.utils.fft(pos, psire[:, k, 1], pmod=pmod)['a0']
           for k in range(np.shape(psire)[1]-1, -1, -1)]
    B0 = make_interp_spline(ire[::-1,1], B0i)
    alfa0_ai = [femagtools.utils.fft(pos, psire[:, k, 0], pmod=pmod)['alfa0']
                for k in range(np.shape(psire)[1])]
    alfa0_a = make_interp_spline(ire[:,0], alfa0_ai)
    alfa0_bi = [femagtools.utils.fft(pos, psire[:, k, 1], pmod=pmod)['alfa0']
                for k in range(np.shape(psire)[1]-1, -1, -1)]
    alfa0_b = make_interp_spline(ire[::-1,1], alfa0_bi)

    Tqi = [femagtools.utils.fft(pos, torq[:, k], pmod=pmod)['a']
           for k in range(np.shape(torq)[1])]
    Tq = make_interp_spline(ire[:, 0], Tqi)
    Tq0i = [femagtools.utils.fft(pos, torq[:, k], pmod=pmod)['a0']
            for k in range(np.shape(torq)[1])]
    Tq0 = make_interp_spline(ire[:, 0], Tq0i)
    alfa0_t = [femagtools.utils.fft(pos, torq[:, k], pmod=pmod)['alfa0']
               for k in range(np.shape(torq)[1])]

    T0 = np.mean([femagtools.utils.fft(pos, psire[:, k, 0], pmod=pmod)['T0']
                  for k in range(np.shape(psire)[1])])
    pp = 360/T0

    def torque(phi, i):
        try:
            alfa0 = np.ones(len(i))*np.mean(alfa0_t)
            alfa0[i < 0] = alfa0_t[0]
            alfa0[i > 0] = alfa0_t[-1]
        except TypeError:
            alfa0 = np.mean(alfa0_t)
            if i < 0:
                alfa0 = alfa0_t[0]
            if i > 0:
                alfa0 = alfa0_t[-1]
        return Tq(i)*np.cos(pp*phi+alfa0) + Tq0(i)

    def psia(phi, i):
        return A(i)*np.cos(pp*phi+alfa0_a(i))+A0(i)

    def psib(phi, i):
        return B(i)*np.cos(pp*phi+alfa0_b(i))+B0(i)

    def dpsiadi(phi,i):
        return A(i, nu=1)*np.cos(pp*phi+alfa0_a(i))+A0(i,nu=1)
    def dpsiadphi(phi,i):
        return -pp*A(i)*np.sin(pp*phi+alfa0_a(i))
    def dpsibdi(phi,i):
        return B(i, nu=1)*np.cos(pp*phi+alfa0_b(i))+B0(i,nu=1)
    def dpsibdphi(phi,i):
        return -pp*B(i)*np.sin(pp*phi+alfa0_b(i))

    speed = simulation['speed']
    r1 = simulation['r1']
    l1s = simulation.get('l1s',0)
    wm = 2*np.pi*speed
    w1 = pp*wm

    def didt(t, y):
        return [((2*r1*y[0] + wm*(
            dpsiadphi(y[1],y[0]) - dpsibdphi(y[1],-y[0])))/
            (-dpsiadi(y[1],y[0]) - dpsibdi(y[1],-y[0]) -2*l1s)),
                wm]
    tmin = simulation.get('tstart', 0)
    tmax = simulation.get('simultime', 0.1)
    nsamples = simulation.get('nsamples', 400)
    t = np.linspace(tmin, tmax, nsamples)

    def func(x):
        return B(0)*np.sin(pp*x+alfa0_b(0)) - A(0)*np.sin(pp*x+alfa0_a(0))
    phi0 = so.fsolve(func, [0])[0]

    Y0 = [0, phi0]
    sol = ig.solve_ivp(didt, (t[0], t[-1]), Y0, dense_output=True)
    ia = sol.sol(t).T[:, 0]
    pv = find_peaks_and_valleys(t, ia)
    iap = pv['tp'], pv['ip']
    iav = pv['tv'], pv['iv']
    iac = pv['tpeaks'][-1], pv['peaks'][-1]

    logger.info("Ia %.1f %.1f %.1f (phi0 %.4f)",
                iap[1], iav[1], iac[1], phi0)

    def func(x):
        y = torque(wm*t+phi0+x, ia)
        pv = find_peaks_and_valleys(t, y)
        return pv['peaks'][-1] + pv['valleys'][-1]

    dphi = so.fsolve(func, [0])[0]
    torque = torque(wm*t+phi0+dphi, ia)
    pv = find_peaks_and_valleys(t, torque)
    tp = pv['tp'], pv['ip']
    tv = pv['tv'], pv['iv']
    tc = pv['tpeaks'][-1], pv['peaks'][-1]
    logger.info("Torque %.1f %.1f %.1f (dphi %.4f)",
                tp[1], tv[1], tc[1], dphi)

    # rotor position at maximum current:
    trot = min(iav[0], iap[0])
    phirot = float(wm*trot + phi0)
    logger.debug("phirot %.1f", phirot)
    #phirot = simulation['phi']

    scData = {
        'ia': ia.tolist(),
        'ib': (-ia).tolist(),
        'ic': np.zeros(ia.shape).tolist(),
        'time': t.tolist(),
        'torque': torque.tolist(),
        'speed': speed,
        'ikd': iac[1],
        'tkd': tc[1],
        'iks': iap[1] if iap[1] > abs(iav[1]) else iav[1],
        'tks': tp[1] if tp[1] > abs(tv[1]) else tv[1]
    }
    scData['peakWindingCurrents'] = [float(scData['iks']),
                                     -float(scData['iks']), 0]
    if simulation.get('sim_demagn', 0):
        i1max = np.abs(scData['iks'])
        scData['demag'] = demag(femag, machine, simulation,
                                i1max, phirot,
                                calculationMode='psi-torq-rem')
    return scData

def demag(femag, machine, simulation, i1max, phirot,
          calculationMode='psi-torq-rem'):
    """demag simulation using psi-torq-rem (or psi-torq-rem-rot)"""
    ik = simulation['ik']
    logger.info("Demagnetization processing: i1max %.1f ik %.1f",
                i1max, ik)

    num_steps = simulation.get('num_demag_cur_steps', 12)
    if ik < i1max:
        i1min = simulation.get('i1min', abs(ik/3))
        n1 = min(num_steps-2,
                 int(num_steps*1.4*(ik-i1min)/(i1max-i1min)))
        n2 = num_steps - n1
        xtab = np.linspace(i1min/abs(ik),
                           1 + 2*(1 - i1min/abs(ik))/(n1 - 1),
                           n1 + 1)
        b = (i1min-abs(ik))/np.log(i1min/abs(ik))
        itab0 = abs(ik) + b*np.log(xtab)
        if n2 > 2:
            x0 = (2*itab0[-1]-itab0[-2])/i1max
            a = np.log(1/x0)/n2
            i0 = x0*i1max
            i1tab = np.hstack(
                (itab0,
                 i0*np.exp(a*np.linspace(0, n2+1, n2+2))))
        else:
            di = itab0[-1]-itab0[-2]
            if itab0[-1]+di < 1.1*i1max:
                i1tab = np.hstack(
                    (itab0,
                     np.arange(itab0[-1]+1.2*di, 1.2*i1max, 1.4*di)))
            else:
                i1tab = itab0
    else:
        i1min = simulation.get('i1min', i1max/3)
        xtab = np.linspace(i1min/abs(ik),
                           1 + 2*(1 - i1min/abs(ik))/(num_steps - 1),
                           num_steps + 1)
        b = (i1min - abs(ik))/np.log(i1min/abs(ik))
        a = abs(ik)
        i1tab = 1.1*(a + np.log(xtab)*b)

    if simulation.get('sc_type', 3) == 3:
        curvec = i1tab
    else:
        curvec = [[-a, a, 0] for a in i1tab]
    simulation.update({
        'calculationMode': calculationMode,
        'phi': phirot,
        'magntemp': simulation['Tmag'],
        'curvec': curvec})
    try:
        del simulation['airgap_induc']
    except KeyError:
        pass
    _ = femag(machine, simulation, fslfile='demag.fsl')

    ptr = np.loadtxt(pathlib.Path(femag.workdir) / "psi-torq-rem.dat")
    i1 = np.concat(([0], np.max(
        np.abs(ptr[:,1:4]), axis=1)))
    rr = np.concat(([1], ptr[:,-1]))
    dmag = {'Hk': simulation['Hk'],
            'Tmag': simulation['Tmag'],
            'i1': i1.tolist(),
            'i1max': float(np.abs(i1max)),
            'rr': rr.tolist()}
    # critical current
    try:
        y = 0.95
        if 0.95 < rr[-1] < 0.97:
            # extrapolate
            x0, x1 = i1[-2], i1[-1]
            y0, y1 = rr[-2], rr[-1]
        else:
            k = np.where(rr < 0.95)[0][0]
            x0, x1 = i1[k-1], i1[k]
            y0, y1 = rr[k-1], rr[k]
        m = (y1-y0)/(x1-x0)
        dmag['i1c'] = (y - y0)/m + x0
    except IndexError:
        pass
    if abs(i1max) < i1[-1]:
        rrf = make_interp_spline(i1, rr)
        dmag['rr_i1max'] = float(rrf(abs(i1max)))
    return dmag
