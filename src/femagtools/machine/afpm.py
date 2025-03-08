""" Axial Flux PM Machine

"""
import logging
import numpy as np
from pathlib import Path
import shutil
from .. import poc
from .. import parstudy
from .. import model
from .. import utils
from .. import windings
from .. import femag
from .. import ecloss
from .. import nc
from .. import bch

from scipy.interpolate import make_interp_spline, RegularGridInterpolator, RectBivariateSpline
from scipy.integrate import quad
import copy
from matplotlib.colors import to_rgb
from multiprocessing import Pool

logger = logging.getLogger(__name__)

AFM_TYPES = (
    "S1R1",      # 1 stator, 1 rotor
    "S1R2",      # 1 stator, 2 rotor, 1 half simulated
    "S1R2_all",  # 1 stator, 2 rotor, all simulated
    "S2R1",      # 2 stator, 1 rotor, 1 half simulated
    "S2R1_all"   # 2 stator, 1 rotor, all simulated
)

def stator_template_name(machine):
    for k in machine['stator']:
        if type(machine['stator'][k]) == dict:
            return k
    return ''

def magnet_template_name(machine):
    for k in machine['magnet']:
        if type(machine['magnet'][k]) == dict:
            return k
    return ''

def jsonify(v):
    vtype = type(v)
    if vtype == np.float32:
        return float(v)
    if vtype == np.ndarray:
        return v.tolist()
    if vtype == dict:
        for k in v:
            v[k] = jsonify(v[k])
    return v


def get_magdata(task):
    basedir = Path(task.directory)
    bchfile_list = sorted(basedir.glob(
        '*_[0-9][0-9][0-9].B*CH'))
    if bchfile_list:
        result = bch.Reader()
        with open(bchfile_list[-1]) as f:
            logger.debug("Reading %s",
                         bchfile_list[-1])
            result.read(f)
    fnc = list(basedir.glob('*.nc'))[0]
    ibeta = len(result.linearForce)-1
    ncmod = nc.read(fnc)
    result.losses[0]['magnet_data'] = [jsonify(pm)
                                       for pm in ncmod.get_magnet_data(
                                               ibeta=ibeta)]
    return result


def num_agnodes(Q, p, pw, ag):
    """return total number of nodes in airgap per pole
     Args:
      Q: (int) number of slots
      p: (int) number of pole pairs
      pw: (float) pole width
      ag: (float) airgap height
    """
    num_nodes = np.arange(6, 120, 6)
    i = np.argmin(np.abs(pw/num_nodes - ag/2))
    nag = num_nodes[i-1]
    if p*nag % Q:
        lcm = np.lcm(Q, 2*p)//p
        nmin, nmax = 1, num_nodes[-1]//lcm
        num_nodes = np.array(
            [i*lcm for i in range(nmin, nmax) if i*lcm % 6 == 0])
        i = np.argmin(np.abs(pw/num_nodes - ag/2))
        if i > 0:
            nag = num_nodes[i-1]
        else:
            nag = num_nodes[i]
        # TODO nodedist 0.5, 2, 4, 6
    return nag


def _integrate(radius, pos, val):
    interp = RegularGridInterpolator((radius, pos), val)
    def func(x, y):
        return interp((x, y))
    return [quad(func, radius[0], radius[-1],
                  args=(p,), limit=200)[0]
            for p in pos]


def _integrate1d(radius, val):
    interp = make_interp_spline(radius, val, k=1)
    def func(x):
        return interp((x))
    return quad(func, radius[0], radius[-1], limit=100)[0]


def ld_interpol(i1, beta, v):
    '''interpolate Ld at beta angle 0°, -180°'''
    # ld
    cur = copy.deepcopy(i1)
    betad = copy.deepcopy(beta)
    if np.amin(beta) < -90 and \
        np.amax(beta) > -90:
        # motor and generator
        v[0] = v[1]
        v[-1] = v[-2]

        dbeta = np.abs(beta[0][0] - beta[1][0])
        bp = [[beta[0][0]-dbeta for i in range(len(np.unique(i1)))]] + beta[1:-1] + \
            [[dbeta for i in range(len(np.unique(i1)))]]
    else:
        v[-1] = v[-2]
        dbeta = np.abs(beta[0][0] - beta[1][0])
        bp = beta[0:-1] + \
            [[dbeta for i in range(len(np.unique(i1)))]]

    return RectBivariateSpline(np.unique(bp), np.unique(cur),
                               np.array(v)).ev(*[betad, i1]).tolist()


def lq_interpol(i1, beta, v):
    '''interpolate Lq at beta -90°'''
    if -90 not in np.unique(beta):
        return v
    # lq
    betad = copy.deepcopy(beta)
    if np.amin(beta) < -90 and \
        np.amax(beta) > -90:
        # motor and generator
        inx = np.argwhere(np.array(beta) == -90).squeeze()
        v.pop(inx[0, 0])
        bp = beta[0:inx[0, 0]] + beta[inx[0, 0]+1:]
        cp = i1[0:inx[0, 0]] + i1[inx[0, 0]+1:]
    else:
        v[0] = v[1]
        dbeta = np.abs(beta[0][0] - beta[1][0])
        bp = [[-90-dbeta for i in i1[0]]] + beta[1::]
        cp = i1
    cur = copy.deepcopy(cp)
    return RectBivariateSpline(np.unique(bp), np.unique(cur), \
         np.array(v)).ev(*[betad, i1]).tolist()


def parident(workdir, engine, temp, machine,
             magnetizingCurves, magnetMat=[], condMat=[],
             **kwargs):
    """return list of parameters of equivalent circuit for AFPM

    arguments:
    workdir: (str) directory for intermediate files
    engine: (object) calculation driver (multiproc, docker, condor)

    temp: list of magnet temperatures in degree Celsius
    machine: dict() with machine parameters
    magnetizingCurves: list of dict() with BH curves
    magnetMat: list of dict() with magnet material properties
    condMat: list of dict() with conductor material properties

    optional arguments:
    num_slices: number of slices (default 3)
    num_cur_steps: number of current steps (default 5)
    num_beta_steps: number of current steps (default 7 per quadrant)
    speed: rotor speed in 1/s (default 160/p)
    i1_max: maximum current in A rms (default approx 3*i1nom)
    cmd: femag executable
    """
    try:
        if machine['afmtype'] not in AFM_TYPES:
            raise ValueError(f"invalid afm type {machine['afmtype']}")
    except KeyError:
        raise ValueError("missing key afmtype")

    di = machine['inner_diam']
    Q1 = machine['stator']['num_slots']
    slotmodel = stator_template_name(machine)
    hs = machine['stator'][slotmodel]['slot_height']
    wdgk = 'windings' if 'windings' in machine else 'winding'
    N = machine[wdgk]['num_wires']
    Jmax = 15  # max current density in A/mm2

    i1_max = kwargs.get(
        'i1_max',
        round(0.28*np.pi*hs*(di+hs)/Q1/N*Jmax*1e5)*10 * \
        machine[wdgk].get('num_par_wdgs', 1))

    p = machine['poles']
    slotmodel = magnet_template_name(machine)
    if np.isscalar(machine['magnet'][slotmodel]['rel_magn_width']):
        num_slices = kwargs.get('num_slices', 3)
        rmagw = num_slices*[machine['magnet'][slotmodel]['rel_magn_width']]
    else:
        rmagw = machine['magnet'][slotmodel]['rel_magn_width']
        num_slices = len(rmagw)

    logger.info("num_slices %d rmagw %s", num_slices, rmagw)
    lfe = get_arm_lengths(machine['outer_diam'],
                          machine['inner_diam'],
                          num_slices)
    pole_width = get_pole_widths(machine['outer_diam'],
                                 machine['inner_diam'],
                                 p,
                                 num_slices)

    speed = kwargs.get('speed', 160/p)
    linspeed = [speed*p*pw for pw in pole_width]

    if "num_agnodes" not in machine:
        for pw in pole_width:
            machine['num_agnodes'] = num_agnodes(Q1, p//2, pw,
                                                 machine['airgap'])
    nlparvardef = {
        "decision_vars": [
            {"values": pole_width,
             "name": "pole_width"},
            {"values": rmagw,
             "name": f"magnet.{slotmodel}.rel_magn_width"},
            {"values": lfe,
             "name": "lfe"},
            {"values": linspeed, "name": "speed"}
        ]
    }
    machine['pole_width'] = np.pi * machine['inner_diam']/machine['poles']
    machine['lfe'] = machine['outer_diam'] - machine['inner_diam']

    beta_min = kwargs.get('beta_min', -90)
    default_steps = round(-beta_min/15)+1
    num_cur_steps = kwargs.get('num_cur_steps', 5)
    num_beta_steps = kwargs.get('num_beta_steps', default_steps)
    parvardef = {
        "decision_vars": [
            #            {"values": sorted(2*temp), "name": "magn_temp"},
            {"steps": num_beta_steps, "bounds": [beta_min, 0],
             "name": "angl_i_up"},
            {"steps": num_cur_steps,
             "bounds": [i1_max/num_cur_steps, i1_max], "name": "current"}
        ]
    }

    ldq = []
    for magtemp in temp:
        nlcalc = dict(
            calculationMode="cogg_calc",
            num_move_steps=60,
            magn_temp=magtemp,
            poc=poc.Poc(999),
            speed=0)
        logging.info("Noload simulation")
        if kwargs.get('use_multiprocessing', True):
            pstudy = parstudy.List(
                    workdir, condMat=condMat, magnets=magnetMat,
                    magnetizingCurves=magnetizingCurves,
                    cmd=kwargs.get('cmd', None))

            nlresults = pstudy(nlparvardef, machine, nlcalc, engine)
            if nlresults['status'].count('C') != len(nlresults['status']):
                raise ValueError(
                    f"Noload simulation failed {nlresults['status']}")
        else:
            nlresults = {"x": [], "f": []}
            i = 0

            for pw, le, sp, rmw in zip(pole_width, lfe, linspeed, rmagw):
                nlmachine = {k: machine[k] for k in machine}
                nlmachine['pole_width'] = pw
                nlmachine['magnet'][slotmodel]['rel_magn_width'] = rmw
                nlmachine['lfe'] = le
                nlcalc.update({"speed": sp})
                nlsubdir = f'{workdir}/{i}'
                nlworkdir = Path(nlsubdir)
                if nlworkdir.exists():
                    shutil.rmtree(nlworkdir)
                nlworkdir.mkdir(exist_ok=True)
                noloadsim = femag.Femag(
                    nlworkdir, condMat=condMat, magnets=magnetMat,
                    magnetizingCurves=magnetizingCurves,
                    cmd=kwargs.get('cmd', None))
                r = noloadsim(nlmachine, nlcalc)
                nlresults['f'].append({k: v for k, v in r.items()})
                i = i + 1
        nlresults.update(process(lfe, pole_width, machine, nlresults['f']))
        weights = nlresults['weights']
        current_angles = nlresults['f'][0]['current_angles']
        results = []
        i = 0
        for l, pw, rmw in zip(lfe, pole_width, rmagw):
            mpart = {k: machine[k] for k in machine if k != slotmodel}
            mpart['pole_width'] = pw
            mpart['lfe'] = l
            mpart['magnet'][slotmodel]['rel_magn_width'] = rmw
            subdir = f"{workdir}/{i}"

            simulation = dict(
                calculationMode="torq_calc",
                wind_temp=20.0,
                magn_temp=magtemp,
                angl_i_up=0.0,
                magn_height=machine['magnet'][slotmodel]['magn_height'],
                yoke_height=machine['magnet'][slotmodel].get(
                    'yoke_height', 0),
                current=1,
                poc=poc.Poc(999,
                            parameters={
                                'phi_voltage_winding': current_angles}),
                num_move_steps=60,
                speed=linspeed[i],
                num_par_wdgs=machine[wdgk].get('num_par_wdgs', 1))

            if kwargs.get('use_multiprocessing', True):
                gpstudy = parstudy.Grid(
                    subdir, condMat=condMat, magnets=magnetMat,
                    magnetizingCurves=magnetizingCurves,
                    cmd=kwargs.get('cmd', None),
                    result_func=get_magdata)
                lresults = gpstudy(parvardef, mpart, simulation, engine)

            else:
                lresults = {"x": [], "f": []}
                domain_beta = np.linspace(beta_min, 0, num_beta_steps).tolist()
                domain_cur = np.linspace(i1_max/num_cur_steps, i1_max, num_cur_steps).tolist()
                dir_index = 0
                for cur in domain_cur:
                    for be in domain_beta:
                        simulation['angl_i_up'] = be
                        simulation['current'] = cur
                        lresults['x'].append([be, cur])
                        subsubdir = subdir + f'/{dir_index}'
                        dir_index = dir_index + 1
                        lworkdir = Path(subsubdir)
                        if lworkdir.exists():
                            shutil.rmtree(lworkdir)
                        lworkdir.mkdir(exist_ok=True)
                        loadsim = femag.Femag(
                            lworkdir, condMat=condMat, magnets=magnetMat,
                            magnetizingCurves=magnetizingCurves,
                            cmd=kwargs.get('cmd', None))
                        r = loadsim(mpart, simulation)
                        lresults['f'].append({k: v for k, v in r.items()})

            f = [{k: bch[k]
                  for k in ('linearForce', 'flux', 'losses', 'lossPar')}
                 for bch in lresults['f']]
            lresults['f'] = f
            results.append(lresults)
            i += 1

        postp = []
        if kwargs.get('use_multiprocessing', True):
            with Pool() as p:
                for r in p.starmap(process,
                                   [(lfe, pole_width, machine, bch)
                                    for bch in zip(*[r['f']
                                                     for r in results])]):
                    torque = np.mean(r.pop('torque'))
                    r['torque'] = torque
                    r.update(_psidq_ldq(r, nlresults))
                    postp.append(r)
        else:
            for r in [process(lfe, pole_width, machine, bch)
                      for bch in zip(*[r['f'] for r in results])]:
                torque = np.mean(r.pop('torque'))
                r['torque'] = torque
                r.update(_psidq_ldq(r, nlresults))
                postp.append(r)

        r1 = postp[0]['r1']
        i1 = [r['i1'] for r in postp][::num_beta_steps]
        beta = [r['beta'] for r in postp][:num_beta_steps]
        # Note: Psi RMS Values aligned with BCH/BATCH
        psid = np.reshape([r['psid'] for r in postp],
                          (-1, num_beta_steps)).T/np.sqrt(2)
        psiq = np.reshape([r['psiq'] for r in postp],
                          (-1, num_beta_steps)).T/np.sqrt(2)

        ld = np.reshape([r['Ld'] for r in postp],
                          (-1, num_beta_steps)).T.tolist()
        lq = np.reshape([r['Lq'] for r in postp],
                          (-1, num_beta_steps)).T.tolist()
        # interpolation ld, lq
        curr, angl = [], []
        for cr in range(len(beta)):
            curr.append(i1)
        for al in beta:
            tmp = []
            for cr in range(len(i1)):
                tmp.append(al)
            angl.append(tmp)
        try:
            xx, yy = copy.deepcopy(curr), copy.deepcopy(angl)
            ld = ld_interpol(xx, yy, ld)
            xx, yy = copy.deepcopy(curr), copy.deepcopy(angl)
            lq = lq_interpol(xx, yy, lq)
        except:
            ld = np.zeros_like(psid).tolist()
            lq = np.zeros_like(psid).tolist()

        torque = np.reshape([r['torque'] for r in postp],
                            (-1, num_beta_steps)).T
        losses = {k: np.flip(np.reshape([r['plfe'][k] for r in postp],
                                        (-1, num_beta_steps)),
                             axis=1).T.tolist()
                  for k in postp[0]['plfe']}
        losses.update({'speed': speed,
                       'hf': postp[0]['lossPar']['hf'],
                       'ef': postp[0]['lossPar']['ef'],
                       'magnet': np.flip(
                           np.reshape([r['plmag'] for r in postp],
                                      (-1, num_beta_steps)),
                           axis=1).T.tolist()})
        ldq.append({'temperature': magtemp,
                    'i1': i1, 'beta': beta,
                    'psid': psid.tolist(), 'psiq': psiq.tolist(),
                    'ld': ld, 'lq': lq,
                    'torque': torque.tolist(),
                    'losses': losses})
        # T = 3/2 p (Psid iq - Psiq id)
        #iq, id = femagtools.machine.utils.iqd(*np.meshgrid(beta, i1))

    return {'m': machine[wdgk]['num_phases'],
            'p': machine['poles']//2, 'weights': weights,
            'rotor_mass': sum(weights[1]), "kfric_b": 1,
            'ls1': 0, 'r1': r1, 'ldq': ldq}


def process(lfe, pole_width, machine, bch):
    """process results: torque, voltage (emf), losses

    Arguments:
    lfe: (float) active machine length
    pole_width: (float) pole width
    machine: (dict) machine
    bch: (list of dict) linearForce, flux, losses

    Returns dict with keys:
    pos: (list of float)
    r1: winding resistance
    torque: (list of float)
    emf: (list of float)
    emf_amp: (float)
    emf_angle: (float)
    freq: (float)
    currents: (list of list of float)
    plfe, plmag, plcu: (list of float)
    """
    model_type = machine['afmtype']
    num_slots = machine['stator']['num_slots']
    mmod = model.MachineModel(machine)
    slots_gen = mmod.stator['num_slots_gen']
    scale_factor = _get_scale_factor(model_type, num_slots, slots_gen)
    endpos = [2*pw*1e3 for pw in pole_width]
    displ = [[d for d in r['linearForce'][-1]['displ']
              if d < e*(1+1/len(r['linearForce'][-1]['displ']))]
             for e, r in zip(endpos, bch)]
    radius = [pw*machine['poles']/2/np.pi for pw in pole_width]
    rotpos = [np.array(d)/r/1000 for r, d in zip(radius, displ)]
    n = len(rotpos[0])
    currents = [bch[0]['flux'][k][0]['current_k'][:n]
                for k in bch[0]['flux']]

    if len(pole_width) > 1:
        # check homogenity:
        if np.diff([len(d) for d in displ]).any():
            raise ValueError(
                f"inhomogenous number of steps: {[len(d) for d in displ]}")
        lfx = [r['linearForce'][-1]['force_x'] for r in bch]
        torque = _integrate(radius, rotpos[0], np.array(
            [r*scale_factor*np.array(fx[:-1])/l
             for l, r, fx in zip(lfe, radius, lfx)]))

        voltage = {k: [scale_factor * np.array(ux[:-1])/l
                       for l, ux in zip(lfe, [r['flux'][k][-1]['voltage_dpsi']
                                              for r in bch])]
                   for k in bch[0]['flux']}
        emf = [_integrate(radius, rotpos[0], np.array(voltage[k]))
               for k in voltage]

        fluxxy = {k: [scale_factor * np.array(flx[:-1])/l
                      for l, flx in zip(lfe, [r['flux'][k][-1]['flux_k']
                                              for r in bch])]
                  for k in bch[0]['flux']}
        flux = [_integrate(radius, rotpos[0], np.array(fluxxy[k]))
                for k in fluxxy]
    else:
        r = radius[0]
        torque = [r*scale_factor*fx
                  for fx in bch[0]['linearForce'][-1]['force_x'][:-1]]
        voltage = {k: [scale_factor * ux
                       for ux in bch[0]['flux'][k][-1]['voltage_dpsi'][:-1]]
                   for k in bch[0]['flux']}
        emf = [voltage[k][:n] for k in voltage]
        fluxxy = {k: [scale_factor * np.array(flx)
                      for flx in bch[0]['flux'][k][-1]['flux_k']]
                  for k in bch[0]['flux']}
        flux = [fluxxy[k][:n] for k in fluxxy]

    pos = (rotpos[0]/np.pi*180)
    emffft = utils.fft(pos, emf[0])

    scale_factor = 1 if model_type == 'S1R1' else 2
    styoke = {}
    stteeth = {'hyst': 0, 'eddy': 0}
    rotor = {}
    for k in ('hyst', 'eddy'):
        if len(pole_width) > 1:
            styoke[k] = _integrate1d(radius, scale_factor*np.array(
                [sum(b['losses'][-1]['stator']['stfe'][k])/l
                 for l, b in zip(lfe, bch)]))
            rotor[k] = _integrate1d(radius, scale_factor*np.array(
                [sum(b['losses'][-1]['rotor']['----'][k])/l
                 for l, b in zip(lfe, bch)]))
        else:
            styoke[k] = scale_factor*sum(
                bch[0]['losses'][-1]['stator']['stfe'][k])
            rotor[k] = scale_factor*sum(
                bch[0]['losses'][-1]['rotor']['----'][k])

    if 'magnet_data' in bch[0]['losses'][0]:
        k = 'magnetH'
        nsegx = machine['magnet'].get('num_segments', 1)
        if type(nsegx) == int:
            nsegx = [nsegx]*len(bch)
        for nx, b in zip(nsegx, bch):
            pm = b['losses'][0]['magnet_data']
            magloss = ecloss.MagnLoss(magnet_data=[pm])
            ml = magloss.calc_losses_ialh2(nsegx=nx)
            b['losses'][-1][k] = ml[0]
    else:
        k = 'magnetJ'
    if len(pole_width) > 1:
        maglosses = _integrate1d(radius, scale_factor*np.array(
            [b['losses'][-1][k]/lz
             for lz, b in zip(lfe, bch)]))
    else:
        maglosses = scale_factor*bch[0]['losses'][-1][k]

    freq = bch[0]['losses'][-1]['stator']['stfe']['freq'][0]

    wdg = windings.Winding(
        dict(Q=mmod.stator['num_slots'],
             p=mmod.poles//2,
             m=mmod.winding['num_phases'],
             l=mmod.winding['num_layers']))
    slotmodel = stator_template_name(machine)
    cufill = mmod.winding.get('cufilfact', 0.4)
    aw = (mmod.stator[slotmodel]['slot_width']*
              mmod.stator[slotmodel]['slot_height']*
              cufill/mmod.winding['num_wires']/mmod.winding['num_layers'])
    r1 = wdg_resistance(wdg, mmod.winding['num_wires'],
                        mmod.winding['num_par_wdgs'],
                        aw,
                        mmod.outer_diam, mmod.inner_diam)
    i1 = np.mean([np.max(c) for c in currents])/np.sqrt(2)
    plcu = mmod.winding['num_phases']*i1**2*r1
    weights = np.array([[0, 0, 0], [0, 0, 0]])
    try:
        for b in bch:
            weights = weights + b['weights']
        weights *= scale_factor
    except KeyError as exc:
        #logger.warning("missing key %s", exc)
        pass
    return {
        'weights': weights.tolist(),
        'pos': pos.tolist(), 'r1': r1,
        'torque': torque,
        'flux_k': flux,
        'emf': emf,
        'emf_amp': emffft['a'], 'emf_angle': emffft['alfa0'],
        'freq': freq,
        'currents': currents,
        'beta': bch[0]['losses'][0]['beta'],
        'plfe': {'styoke_hyst': styoke['hyst'], 'styoke_eddy': styoke['eddy'],
                 'stteeth_hyst': stteeth['hyst'], 'stteeth_eddy': stteeth['eddy'],
                 'rotor_hyst': rotor['hyst'], 'rotor_eddy': rotor['eddy']},
        'plmag': maglosses, 'lossPar': {k: bch[0]['lossPar'][k]
                                        for k in bch[0]['lossPar']},
        'plcu': plcu}


def _psidq_ldq(results, nlresults):
    """ calculate currents, winding fluxes psi and inductances L
    Arguments:
    results: (dict) emf, currents, freq of load simulation
    nlresults: (dict) emf of no load simulation

    returns dict of currents, psi, L
    psid, psiq, psim, i1, id, iq, Ld, Lq
    """
    gamma = -(results['emf_angle'] - nlresults['emf_angle'])
    w1 = 2*np.pi*results['freq']
    psid = np.cos(gamma)*results['emf_amp']/w1
    psiq = -np.sin(gamma)*results['emf_amp']/w1
    psim = nlresults['emf_amp']/w1
    i1 = np.mean([np.max(c)
                  for c in results['currents']])/np.sqrt(2)
    beta = results['beta']/180*np.pi
    id = np.sqrt(2)*i1*np.sin(beta)
    iq = np.sqrt(2)*i1*np.cos(beta)
    Ld = Lq = float('nan')
    if np.abs(id) > 0:
        Ld = (psid - psim)/id
    if np.abs(iq) > 0:
        Lq = psiq/iq
    return {
        'psid': psid,
        'psiq': psiq,
        'psim': psim,
        'i1': i1,
        'beta': results['beta'],
        'id': id,
        'iq': iq,
        'Ld': Ld,
        'Lq': Lq}

def _get_scale_factor(model_type, num_slots, slots_gen):
    """Determines the scale factor
    Arguments:
    model_type : (str) type of model
    num_slots: (int) total number of stator slots
    slots_gen: (int) number of slots in model

    Returns scale factor based on number of poles, number of poles simulated
        and the model type
    """
    segments = num_slots / slots_gen
    if model_type in ("S2R1", "S1R2"):
        return 2 * segments
    return segments


def get_arm_lengths(outer_diam, inner_diam, num_slices):
    d = outer_diam - inner_diam
    if num_slices > 2:
        return [d/(4*(num_slices-1))] + [
            d/(2*(num_slices-1))
            for i in range(1,num_slices-1)] + [d/(4*(num_slices-1))]
    return [d/2]


def get_pole_widths(outer_diam, inner_diam, poles, num_slices):
    d = outer_diam - inner_diam
    if num_slices > 2:
        return [np.pi * inner_diam/poles] + [
            np.pi * (inner_diam + d*i/(num_slices - 1))/poles
            for i in range(1, num_slices-1)] + [np.pi * outer_diam/poles]
    return [np.pi * (outer_diam+inner_diam)/2/poles]


def wdg_resistance(wdg, n, g, aw, outer_diam, inner_diam,
                   sigma=56e6):
    """return winding resistance per phase in Ohm
    Arguments:
    wdg: (Winding) winding
    n: (int) number of wires per coil side
    g: (int) number of parallel coil groups
    aw: wire cross section area m2
    outer_diam, inner_diam: diameters m
    sigma: (float) conductivity of wire material 1/Ohm m
    """
    # mean length of one turn
    lt = 2.4*((outer_diam-inner_diam)+np.pi/wdg.Q*(outer_diam+inner_diam) + 16e-3)
    return wdg.turns_per_phase(n, g)*lt/sigma/aw/g


def _get_copper_losses(scale_factor, bch):
    """return copper losses from bch files"""
    try:
        wdgk = 'winding'
        cu_losses = sum([b['losses'][0][wdgk] for b in bch])
        return scale_factor*cu_losses
    except KeyError:
        return 0  # noload calc has no winding losses


def _set_plot_attributes(ax):
    ax.set_aspect('equal')
    for loc, spine in ax.spines.items():
        spine.set_color('none')  # don't draw spine
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks([])


def _get_colors(colors, delta):
    if delta == 0.0:
        return colors
    new_colors = []
    for col in colors:
        rgb = to_rgb(col)
        r, g, b = rgb
        col = (max(0.0, min(r+delta, 1.0)),
               max(0.0, min(g+delta, 1.0)),
               max(0.0, min(b+delta, 1.0)))
        new_colors.append(col)
    return new_colors


def _draw_vertical_magnets(ax,
                           poles,
                           xr, yr,
                           Rr,
                           xoff, yoff,
                           delta=0.0):
    color = ['green', 'red']
    color = _get_colors(color, delta)
    for i in range(poles):
        ax.fill(xr+xoff, yr+yoff,
                facecolor=color[i%2], edgecolor=color[i%2])
        xr, yr = np.dot(Rr, [xr, yr])
    return


def _draw_vertical_slots(ax,
                         Q,
                         r,
                         alpha,
                         xoff, yoff,
                         delta=0.0):
    color = ['skyblue', 'blue']
    color = _get_colors(color, delta)
    taus = 2*np.pi/Q
    for n in range(Q):
        beta0 = np.linspace(n*taus, n*taus + taus/2-alpha[0], 5)
        beta1 = np.linspace(n*taus, n*taus + taus/2-alpha[1], 5)
        xr = np.concatenate((
            r[0]*np.cos(beta0), r[1]*np.cos(beta1[::-1])))+xoff
        yr = np.concatenate((
            r[0]*np.sin(beta0), r[1]*np.sin(beta1[::-1])))+yoff
        ax.fill(xr, yr, color=color[0])
        beta0 = np.linspace(n*taus + taus/2+alpha[0], (n+1)*taus, 5)
        beta1 = np.linspace(n*taus + taus/2+alpha[1], (n+1)*taus, 5)
        xr = np.concatenate((
            r[0]*np.cos(beta0), r[1]*np.cos(beta1[::-1])))+xoff
        yr = np.concatenate((
            r[0]*np.sin(beta0), r[1]*np.sin(beta1[::-1])))+yoff
        ax.fill(xr, yr, color=color[0])


def vertical_plot(machine, ax):
    """plots afpm stator and rotor (vertical section)
    Args:
      dy1, dy1: float outer, inner diameter
      rel_magn_width: float rel magnet width 0..1
      Q: number of stator slots
      poles: number of poles
      slot_width: width of stator slot
    """
    logger.debug("begin of vertical_plot()")

    model_type = machine['afmtype'][0:4]
    dy1 = machine['outer_diam']*1e3
    dy2 = machine['inner_diam']*1e3
    try:
        rel_magn_width = max(machine['magnet']['afm_rotor']['rel_magn_width'])
    except TypeError:
        rel_magn_width = machine['magnet']['afm_rotor']['rel_magn_width']
    Q = machine['stator']['num_slots']
    slot_width = machine['stator']['afm_stator']['slot_width']*1e3
    poles = machine['poles']

    # prepare Magnets
    theta = np.linspace(np.pi/poles*(1-rel_magn_width),
                        np.pi/poles*(1+rel_magn_width),
                        10)
    xr = np.concatenate((dy1/2 * np.cos(theta), dy2/2 * np.cos(theta[::-1])))
    yr = np.concatenate((dy1/2 * np.sin(theta), dy2/2 * np.sin(theta[::-1])))
    rtheta = 2*np.pi/poles
    Rr = np.array([
        [np.cos(rtheta), -np.sin(rtheta)],
        [np.sin(rtheta),  np.cos(rtheta)]
    ])

    # prepare Slots
    taus = 2*np.pi/Q
    r = np.array([dy2/2, dy1/2])
    alpha = np.arctan2(slot_width/2, r)

    yoff = 0.0
    xoff = 0.0
    y_shift = -2
    x_shift = dy1-dy2
    ma_delta = 0.0  # color
    sl_delta = 0.0  # color

    # Draw
    if model_type in ("S1R2"):  # 2 rotor
        _draw_vertical_magnets(ax, poles, xr, yr, Rr, xoff, yoff, delta=-0.1)
        yoff += y_shift
        xoff += x_shift

    if model_type in ("S2R1"):  # 2 stator
        sl_delta = -0.1

    _draw_vertical_slots(ax, Q, r, alpha, xoff, yoff, delta=sl_delta)
    yoff += y_shift
    xoff += x_shift

    if model_type in ("S1R2"):  # 2 rotor
        ma_delta = 0.1

    _draw_vertical_magnets(ax, poles, xr, yr, Rr, xoff, yoff, delta=ma_delta)
    yoff += y_shift
    xoff += x_shift

    if model_type in ("S2R1"):  # 2 stator
        sl_delta = 0.0
        _draw_vertical_slots(ax, Q, r, alpha, xoff, yoff, delta=sl_delta)

    _set_plot_attributes(ax)
    logger.debug("end of vertical_plot()")


IRON_NO = 0
IRON_UP = 1
IRON_DOWN = 2

def _draw_horizontal_magnets(ax,
                             poles,
                             magn_height,
                             magn_width,
                             yoke_height,
                             Q,
                             g,
                             taus,
                             dy2,
                             yoff=0.0,
                             iron=IRON_NO
                             ):
    color = ['green', 'red']
    xy = (0, Q//g*taus, Q//g*taus, 0)

    if iron == IRON_UP:
        yy = (yoff-yoke_height,
              yoff-yoke_height,
              yoff,
              yoff)
        yoff -= yoke_height
        ax.fill(xy, yy, color='skyblue')

    taum = dy2*np.pi/poles
    ym = np.array([yoff-magn_height,
                   yoff-magn_height,
                   yoff,
                   yoff])
    yoff -= magn_height

    for n in range(poles//g):
        xl = taum*n + taum*(1 - magn_width)
        xr = taum*(n + 1) - taum*(1 - magn_width)
        xm = (xl, xr, xr, xl)
        ax.fill(xm, ym, color=color[n%2])

    if iron == IRON_DOWN:
        yy = (yoff-yoke_height,
              yoff-yoke_height,
              yoff,
              yoff)
        yoff -= yoke_height
        ax.fill(xy, yy, color='skyblue')
    return yoff


TOOTH_UP = 0
TOOTH_DOWN = 1
TOOTH_ONLY = 2


def _draw_horizontal_slots(ax,
                           slot_height, slot_width, yoke_height,
                           Q, g, taus,
                           yoff=0.0,
                           tooth=TOOTH_DOWN):
    if not tooth == TOOTH_ONLY:
        xx = (0, Q//g*taus, Q//g*taus, 0)
        if tooth == TOOTH_DOWN:
            yy = (yoff-yoke_height,
                  yoff-yoke_height,
                  yoff,
                  yoff)
        else:
            yy = (yoff-slot_height,
                  yoff-slot_height,
                  yoff-slot_height-yoke_height,
                  yoff-slot_height-yoke_height)
        ax.fill(xx, yy, color='skyblue')

    yt = (yoff-slot_height-yoke_height,
          yoff-slot_height-yoke_height,
          yoff, yoff)
    for n in range(Q//g):
        xt = np.array((n*taus, n*taus+(taus-slot_width)/2,
                      n*taus+(taus-slot_width)/2, n*taus))
        ax.fill(xt, yt, color='skyblue')
        xt += slot_width + (taus-slot_width)/2
        ax.fill(xt, yt, color='skyblue')
    return yoff - slot_height - yoke_height


def horizontal_plot(machine, ax):
    logger.debug("begin of horizontal_plot()")

    model_type = machine['afmtype'][0:4]
    dy1 = machine['outer_diam']*1e3
    dy2 = machine['inner_diam']*1e3
    try:
        rel_magn_width = max(machine['magnet']['afm_rotor']['rel_magn_width'])
    except TypeError:
        rel_magn_width = machine['magnet']['afm_rotor']['rel_magn_width']
    magn_height = machine['magnet']['afm_rotor']['magn_height']*1e3
    magn_yoke_height = machine['magnet']['afm_rotor']['yoke_height']*1e3

    Q = machine['stator']['num_slots']
    slot_width = machine['stator']['afm_stator']['slot_width']*1e3
    poles = machine['poles']
    m = 3
    slot_height = machine['stator']['afm_stator']['slot_height']*1e3
    if model_type in ('S2R1', 'S2R1_all'):
        slot_height /= 2
    yoke_height = machine['stator']['afm_stator']['yoke_height']*1e3
    ag = machine['airgap']*1e3

    g = np.gcd(Q, m*poles)//m
    taus = dy2*np.pi/Q

    yoff = 0.0
    if model_type in ('S1R2', 'S1R2_all'):  # 2 rotor
        yoff = _draw_horizontal_magnets(ax, poles,
                                        magn_height, rel_magn_width,
                                        magn_yoke_height,
                                        Q, g, taus, dy2,
                                        yoff=yoff,
                                        iron=IRON_UP)
        yoff -= ag

    tooth = TOOTH_ONLY if model_type in ('S1R2', 'S1R2_all') else TOOTH_DOWN
    yoff = _draw_horizontal_slots(ax,
                                  slot_height, slot_width, yoke_height,
                                  Q, g, taus,
                                  yoff=yoff,
                                  tooth=tooth)
    yoff -= ag

    iron = IRON_DOWN if model_type in ('S1R1', 'S1R2', 'S1R2_all') else IRON_NO
    yoff = _draw_horizontal_magnets(ax, poles,
                                    magn_height, rel_magn_width,
                                    magn_yoke_height,
                                    Q, g, taus, dy2,
                                    yoff=yoff,
                                    iron=iron)
    yoff -= ag

    if model_type in ('S2R1', 'S2R1_all'):  # 2 rotor
        yoff = _draw_horizontal_slots(ax,
                                      slot_height, slot_width, yoke_height,
                                      Q, g, taus,
                                      yoff=yoff,
                                      tooth=TOOTH_UP)

    _set_plot_attributes(ax)
    logger.debug("end of horizontal_plot()")


class AFPM:
    """Axial Flux PM
    Arguments:
    workdir: (str) working directory
    magnetizingCurves: (str) or (object)
    magnetMat: magnet material
    condMat: conductor material
    """
    def __init__(self, workdir, magnetizingCurves='.', magnetMat='',
                 condMat='', cmd=None):
        self.parstudy = parstudy.List(
            workdir, condMat=condMat, magnets=magnetMat,
            magnetizingCurves=magnetizingCurves, cmd=cmd)

    def __call__(self, engine, machine, simulation, num_slices):
        """ run FE simulation

        returns
        pos: (list of float)
        r1: radius
        torque: (list of float)
        emf: (list of float)
        emf_amp: (float)
        emf_angle: (float)
        freq: (float)
        currents: (list of list of float)
        plfe, plmag, plcu: (list of float)
        psid, psiq, psim: (float)
        i1, beta: (float)
        id, iq: (float)
        Ld, Lq: (float)
        """
        try:
            if machine['afmtype'] not in AFM_TYPES:
                raise ValueError(f"invalid afm type {machine['afmtype']}")
        except KeyError:
            raise ValueError("missing key afmtype")
        slotmodel = magnet_template_name(machine)
        logger.info("num_slices %d rmagw %s", num_slices,
                    machine['magnet'][slotmodel]['rel_magn_width'])
        if np.isscalar(machine['magnet'][slotmodel]['rel_magn_width']):
            rmagw = num_slices*[machine['magnet'][slotmodel]['rel_magn_width']]
        else:
            rmagw = machine['magnet'][slotmodel]['rel_magn_width']
            num_slices = len(rmagw)

        lfe = get_arm_lengths(machine['outer_diam'],
                              machine['inner_diam'],
                              num_slices)
        pole_width = get_pole_widths(machine['outer_diam'],
                                     machine['inner_diam'],
                                     machine['poles'],
                                     num_slices)
        linspeed = [simulation['speed']*machine['poles']*pw
                    for pw in pole_width]

        if "num_agnodes" not in machine:
            Q1 = machine['stator']['num_slots']
            p = machine['poles']
            for pw in pole_width:
                machine['num_agnodes'] = num_agnodes(Q1, p//2, pw,
                                                     machine['airgap'])
            logger.info("Num agnodes/pole %d", machine['num_agnodes'])
        parvardef = {
            "decision_vars": [
                {"values": pole_width,
                 "name": "pole_width"},
                {"values": lfe,
                 "name": "lfe"},
                {"values": rmagw,
                 "name": f"magnet.{slotmodel}.rel_magn_width"},
                {"values": linspeed, "name": "speed"}
            ]
        }

        machine['pole_width'] = np.pi * machine['inner_diam']/machine['poles']
        machine['lfe'] = machine['outer_diam'] - machine['inner_diam']
        try:
            machine['magnet'][slotmodel]['rel_magn_width'] = max(
                machine['magnet'][slotmodel]['rel_magn_width'])
        except TypeError:
            pass

        simulation['skew_displ'] = (simulation.get('skew_angle', 0)/180 * np.pi
                                    * machine['inner_diam'])
        nlresults = {}
        if (simulation['calculationMode'] == 'torq_calc'
            and 'poc' not in simulation):
            nlcalc = dict(
                calculationMode="cogg_calc",
                magn_temp=simulation.get('magn_temp', 20),
                num_move_steps=60,
                skew_linear=0,
                skew_steps=0,
                poc=poc.Poc(machine['pole_width']),
                speed=0)
            logging.info("Noload simulation")
            nlresults = self.parstudy(parvardef,
                                      machine, nlcalc, engine)
            if nlresults['status'].count('C') != len(nlresults['status']):
                raise ValueError(f'Noload simulation failed {nlresults["status"]}')
            nlresults.update(process(lfe, pole_width, machine, nlresults['f']))

            current_angles = nlresults['f'][0]['current_angles']
            simulation['poc'] = poc.Poc(
                999,
                parameters={
                    'phi_voltage_winding': current_angles})
            logger.info("Current angles: %s", current_angles)
            simulation['calc_noload'] = 0

        if 'poc' not in simulation:
            simulation['poc'] = poc.Poc(machine['pole_width'])

        self.parstudy.result_func = get_magdata
        logger.info("Simulation %s", simulation['calculationMode'])
        lresults = self.parstudy(
            parvardef, machine, simulation, engine)
        if lresults['status'].count('C') != len(lresults['status']):
            raise ValueError(
                f'{simulation["calculationMode"]} failed {lresults["status"]}')

        results = process(lfe, pole_width, machine, lresults['f'])
        if nlresults:
            results.update(_psidq_ldq(results, nlresults))
        return results
