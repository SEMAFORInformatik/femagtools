import pathlib
import logging
import json
import numpy as np
import scipy.optimize as so
from scipy.interpolate import make_interp_spline
import scipy.integrate as ig
import femagtools.parstudy

logger = logging.getLogger('shortcircuit')

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
        cmd=cmd, result_func=result_func)

def get_shortCircuit_parameters(bch, nload):
    """extracts shortciruit parameters from bch"""
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
            ld=ld,
            lq=lq,
            Hk=bch.magnet['demag_hx'],
            Tmag=bch.magnet['Tmag'],
            psim=psim,
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
    scdata = {}
    calcmode = simulation.get('calculationMode', '')
    simulation.update(
        get_shortCircuit_parameters(bch,
                                    simulation.get('initial', 2)))
    if 'speed' not in simulation:
        simulation['speed'] = bch.dqPar['speed']
    if simulation.get('sc_type', 3) == 3:
        logger.info("3phase short circuit simulation")
        builder = femagtools.fsl.Builder(femag.templatedirs)
        fslcmds = (builder.open_model(femag.model) +
                   builder.create_shortcircuit(simulation))
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
                d = {'displ': [d['displ']
                               for d in bchsc.demag if 'displ' in d],
                     'H_max': [d['H_max']
                               for d in bchsc.demag if 'H_max' in d],
                     'H_av': [d['H_av']
                              for d in bchsc.demag if 'H_av' in d]}
                simulation['i1max'] = bchsc.scData['iks']
                bchsc.scData['demag'] = demag(
                    femag, machine, simulation, engine)
                bchsc.scData['demag'].update(d)
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
        if 'i1max' not in simulation:
            # just a wild guess
            simulation['i1max'] = 4.5*bch.machine['i1']
        logger.info("2phase short circuit simulation i1max = %.0f",
                    simulation['i1max'])
        scdata = shortcircuit_2phase(femag, machine, simulation, engine)

    else:
        logger.warning("Empty shortcircuit results for type %d",
                       simulation.get('sc_type', 'unknown'))
    # must reset calcmode
    if calcmode:
        simulation['calculationMode'] = calcmode
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
    return {'pos': pos.tolist(), 'ire': ire.tolist(),
            'psire': psire.tolist(), 'torq': torq.tolist()}


def shortcircuit_2phase(femag, machine, simulation, engine=0):
    i1max = simulation['i1max']
    num_cur_steps = 4
    i1 = np.linspace(0, i1max, num_cur_steps)
    i1vec = np.concat((-i1[::-1], i1[1:]))
    num_par_wdgs = machine['winding'].get('num_par_wdgs', 1)
    flux_sim = {
        'calculationMode': 'psi-torq-rot',
        'i1max': i1max,
        'curvec': [],
        'num_par_wdgs': num_par_wdgs}

    if engine:
        parstudy = _parstudy_list(femag, sc_result_func)
        parvardef = {
            "decision_vars": [
                {"values": i1vec, "name": "curvec"}]
        }
        results = parstudy(parvardef, machine, flux_sim, engine)

        ire = np.array([r['ire'][0] for r in results['f']])
        pos = np.array(results['f'][0]['pos'])
        phi = pos*np.pi/180
        torq = np.hstack([r['torq'] for r in results['f']])
        psire = np.hstack([r['psire'] for r in results['f']])
    else:
        simulation.update(flux_sim)
        simulation['curvec'] = i1vec.tolist()
        results = femag(machine, simulation)
        class Task:
            def __init__(self, workdir):
                self.directory = workdir
        results = sc_result_func(Task(femag.workdir))
        ire = np.array(results['ire'])
        pos = np.array(results['pos'])
        torq = np.array(results['torq'])
        psire = np.array(results['psire'])

    #with open('results.json', 'w') as fp:
    #    json.dump({'ire': ire.tolist(), 'pos': pos.tolist(),
    #               'torq': torq.tolist(), 'psire': psire.tolist()}, fp)
    logger.info("move steps %d currents %s", len(pos), ire[:,0])

    Ai = [femagtools.utils.fft(pos, psire[:, k, 0])['a']
          for k in range(np.shape(psire)[1])]
    A = make_interp_spline(ire[:,0], Ai)
    A0i = [femagtools.utils.fft(pos, psire[:, k, 0])['a0']
           for k in range(np.shape(psire)[1])]
    A0 = make_interp_spline(ire[:,0], A0i)
    Bi = [femagtools.utils.fft(pos, psire[:, k, 1])['a']
          for k in range(np.shape(psire)[1]-1, -1, -1)]
    B = make_interp_spline(ire[::-1,1], Bi)
    B0i = [femagtools.utils.fft(pos, psire[:, k, 1])['a0']
           for k in range(np.shape(psire)[1]-1, -1, -1)]
    B0 = make_interp_spline(ire[::-1,1], B0i)
    alfa0_ai = [femagtools.utils.fft(pos, psire[:, k, 0])['alfa0']
                for k in range(np.shape(psire)[1])]
    alfa0_a = make_interp_spline(ire[:,0], alfa0_ai)
    alfa0_bi = [femagtools.utils.fft(pos, psire[:, k, 1])['alfa0']
                for k in range(np.shape(psire)[1]-1, -1, -1)]
    alfa0_b = make_interp_spline(ire[::-1,1], alfa0_bi)

    Tqi = [femagtools.utils.fft(pos, torq[:, k])['a']
           for k in range(np.shape(torq)[1])]
    Tq = make_interp_spline(ire[:, 0], Tqi)
    Tq0i = [femagtools.utils.fft(pos, torq[:, k])['a0']
            for k in range(np.shape(torq)[1])]
    Tq0 = make_interp_spline(ire[:, 0], Tq0i)
    alfa0_t = [femagtools.utils.fft(pos, torq[:, k])['alfa0']
               for k in range(np.shape(torq)[1])]

    T0 = np.mean([femagtools.utils.fft(pos, psire[:, k, 0])['T0']
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
    scData['peakWindingCurrents'] = [scData['iks'],
                                     -scData['iks'], 0]
    if simulation.get('sim_demagn', 0):
        scData['demag'] = demag(femag, machine, simulation, engine)
    return scData

def dm_result_func(task):
    basedir = pathlib.Path(task.directory)
    i1rr = []
    for f in sorted(basedir.glob('psi-torq-rem-rot-*.dat')):
        ptr = np.loadtxt(f)
        i1rr.append((np.max(ptr.T[1:4]), np.min(ptr.T[-1])))
    return i1rr

def demag(femag, machine, simulation, engine=0):
    """demag simulation using psi-torq-rem-rot"""
    logger.info("Demagnetization processing")
    i1max = simulation['i1max']
    i1min = simulation.get('i1min', i1max/4)
    num_steps = 7
    b = (i1min-i1max)/np.log(i1min/i1max)
    a = i1max/b
    i1tab = [b*(a+np.log(x))
             for x in np.linspace(i1min/i1max, 1,
                                  num_steps)]

    if simulation.get('sc_type', 3) == 3:
        curvec = [[-a/2, a, -a/2] for a in i1tab]
    else:
        curvec = [[a, -a, 0] for a in i1tab]
    simulation.update({
        'calculationMode': 'psi-torq-rem-rot',
        'curvec': curvec})
    if engine:
        parstudy = _parstudy_list(femag, dm_result_func)
        parvardef = {
            "decision_vars": [
                {"values": curvec, "name": "curvec"}]
        }
        results = parstudy(parvardef, machine, simulation, engine)
        i1rr = np.vstack(
            ((0, 1),
             np.array(results['f']).reshape((-1, 2))))
    else:
        class Task:
            def __init__(self, workdir):
                self.directory = workdir
        _ = femag(machine, simulation)
        i1rr = np.vstack(
            [(0, 1), dm_result_func(Task(femag.workdir))])
    i1, rr = np.array(i1rr).T
    dmag = {'Hk': simulation['Hk'],
            'Tmag': simulation['Tmag'],
            'i1': i1.tolist(),
            'rr': rr.tolist()}
    # critical current
    if np.min(rr) < 0.99:
        k = np.where(rr < 0.99)[0][0]
        dmag['i1c'] = i1[k]
    return dmag
