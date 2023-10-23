""":mod:`femagtools.machine.afpm` -- Axial Flux PM Machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""
import logging
import numpy as np
from pathlib import Path
from .. import poc
from .. import parstudy
from .. import model
from .. import utils
from .. import windings
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.integrate import quad

logger = logging.getLogger(__name__)

AFM_TYPES = (
    "S1R1",      # 1 stator, 1 rotor
    "S1R2",      # 1 stator, 2 rotor, 1 half simulated
    "S1R2_all",  # 1 stator, 2 rotor, all simulated
    "S2R1",      # 2 stator, 1 rotor, 1 half simulated
    "S2R1_all"   # 2 stator, 1 rotor, all simulated
)

def _integrate(radius, pos, val):
    interp = RegularGridInterpolator((radius, pos), val)
    def func(x, y):
        return interp((x, y))
    return [quad(func, radius[0], radius[-1], args=(p,))[0]
            for p in pos]


def _integrate1d(radius, val):
    interp = interp1d(radius, val)
    def func(x):
        return interp((x))
    return quad(func, radius[0], radius[-1])[0]


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
    slotmodel = 'afm_stator'
    hs = machine['stator'][slotmodel]['slot_height']
    wdgk = 'windings' if 'windings' in machine else 'winding'
    N = machine[wdgk]['num_wires']
    Jmax = 15  # max current density in A/mm2

    i1_max = kwargs.get(
        'i1_max',
        round(0.28*np.pi*hs*(di+hs)/Q1/N*Jmax*1e5)*10 * \
        machine[wdgk].get('num_par_wdgs', 1))

    p = machine['poles']
    num_slices = kwargs.get('num_slices', 3)
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
            machine['num_agnodes'] = 6*round(pw/machine['airgap']/4)

            #if machine['num_agnodes'] < nper:
            #    machine['num_agnodes'] = 8*round(pw/machine['airgap']/4)

    nlparvardef = {
        "decision_vars": [
            {"values": pole_width,
             "name": "pole_width"},
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
            {"steps": num_beta_steps, "bounds": [beta_min, 0], "name": "angl_i_up"},
            {"steps": num_cur_steps,
             "bounds": [i1_max/num_cur_steps, i1_max], "name": "current"}
        ]
    }

    pstudy = parstudy.List(
        workdir, condMat=condMat, magnets=magnetMat,
        magnetizingCurves=magnetizingCurves,
        cmd=kwargs.get('cmd', None))

    ldq = []
    for magtemp in temp:
        nlcalc = dict(
            calculationMode="cogg_calc",
            num_move_steps=60,
            magn_temp=magtemp,
            speed=0)
        logging.info("Noload simulation")
        nlresults = pstudy(nlparvardef, machine, nlcalc, engine)
        nlresults.update(process(lfe, pole_width, machine, nlresults['f']))
        current_angles = nlresults['f'][0]['current_angles']

        results = []
        i = 0
        for l, pw in zip(lfe, pole_width):
            mpart = {k: machine[k] for k in machine if k != 'afm_rotor'}
            mpart['pole_width'] = pw
            mpart['lfe'] = l
            subdir = f"{workdir}/{i}"

            gpstudy = parstudy.Grid(
                subdir, condMat=condMat, magnets=magnetMat,
                magnetizingCurves=magnetizingCurves,
                cmd=kwargs.get('cmd', None))

            simulation = dict(
                calculationMode="torq_calc",
                wind_temp=20.0,
                magn_temp=magtemp,
                angl_i_up=0.0,
                magnet_loss=True,
                magn_height=machine['magnet']['afm_rotor']['magn_height'],
                yoke_height=machine['magnet']['afm_rotor'].get(
                    'yoke_height', 0),
                current=1,
                poc=poc.Poc(999,
                            parameters={
                                'phi_voltage_winding': current_angles}),
                num_move_steps=60,
                speed=linspeed[i],
                num_par_wdgs=machine[wdgk].get('num_par_wdgs', 1))

            lresults = gpstudy(parvardef, mpart, simulation, engine)
            f = [{k: bch[k]
                  for k in ('linearForce', 'flux', 'losses', 'lossPar')}
                 for bch in lresults['f']]
            lresults['f'] = f
            results.append(lresults)
            i += 1

        postp = []
        for results in [process(lfe, pole_width, machine, bch)
                        for bch in zip(*[r['f'] for r in results])]:
            torque = np.mean(results.pop('torque'))
            results['torque'] = torque
            results.update(_psidq_ldq(results, nlresults))
            postp.append(results)

        r1 = postp[0]['r1']
        i1 = [r['i1'] for r in postp][::num_beta_steps]
        beta = [r['beta'] for r in postp][:num_beta_steps]
        # Note: Psi RMS Values aligned with BCH/BATCH
        psid = np.reshape([r['psid'] for r in postp],
                          (-1, num_beta_steps)).T/np.sqrt(2)
        psiq = np.reshape([r['psiq'] for r in postp],
                          (-1, num_beta_steps)).T/np.sqrt(2)
        torque = np.reshape([r['torque'] for r in postp],
                            (-1, num_beta_steps)).T
        losses = {k: np.flip(np.reshape([r['plfe'][k] for r in postp],
                                        (-1, num_beta_steps)),
                             axis=1).T.tolist()
                  for k in postp[0]['plfe']}
        losses.update({'speed': speed,
                       'hf': postp[0]['lossPar']['hf'],
                       'ef': postp[0]['lossPar']['ef'],
                       'magnet':np.flip(
                           np.reshape([r['plmag'] for r in postp],
                                      (-1, num_beta_steps)),
                           axis=1).T.tolist()})
        ldq.append({'temperature': magtemp,
                    'i1':i1, 'beta':beta,
                    'psid': psid.tolist(), 'psiq': psiq.tolist(),
                    'torque': torque.tolist(),
                    'losses': losses})
        # T = 3/2 p (Psid iq - Psiq id)
        #iq, id = femagtools.machine.utils.iqd(*np.meshgrid(beta, i1))

    return {'m': machine[wdgk]['num_phases'],
            'p': machine['poles']//2,
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
    scale_factor =_get_scale_factor(model_type, num_slots, slots_gen)
    endpos = [2*pw*1e3 for pw in pole_width]
    displ = [[d for d in r['linearForce'][0]['displ']
              if d < e*(1+1/len(r['linearForce'][0]['displ']))]
             for e, r in zip(endpos, bch)]
    radius = [pw*machine['poles']/2/np.pi for pw in pole_width]
    rotpos = [np.array(d)/r/1000 for r, d in zip(radius, displ)]
    n = len(rotpos[0])
    currents = [bch[0]['flux'][k][0]['current_k'][:n]
                for k in bch[0]['flux']]
    if len(pole_width) > 1:
        torque = _integrate(radius, rotpos[0], np.array(
            [r*scale_factor*np.array(fx[:-1])/l
             for l, r, fx in zip(lfe, radius,
                                 [r['linearForce'][0]['force_x']
                                  for r in bch])]))

        voltage = {k: [scale_factor * np.array(ux[:-1])/l
                       for l, ux in zip(lfe, [r['flux'][k][0]['voltage_dpsi']
                                              for r in bch])]
                   for k in bch[0]['flux']}
        emf = [_integrate(radius, rotpos[0], np.array(voltage[k]))
               for k in voltage]
    else:
        r = radius[0]
        torque = [r*scale_factor*fx
                  for fx in bch[0]['linearForce'][0]['force_x'][:-1]]
        voltage = {k: [scale_factor * ux
                       for ux in bch[0]['flux'][k][0]['voltage_dpsi'][:-1]]
                   for k in bch[0]['flux']}
        emf = [voltage[k][:n] for k in voltage]

    pos = (rotpos[0]/np.pi*180)
    emffft = utils.fft(pos, emf[0])

    styoke = {}
    stteeth = {'hyst':0, 'eddy':0}
    rotor = {}
    for k in ('hyst', 'eddy'):
        if len(pole_width) > 1:
            styoke[k] = _integrate1d(radius, scale_factor*np.array(
                [sum(b['losses'][0]['stator']['stfe'][k])/l
                 for l, b in zip(lfe, bch)]))
            rotor[k] = _integrate1d(radius, scale_factor*np.array(
                [sum(b['losses'][0]['rotor']['----'][k])/l
                 for l, b in zip(lfe, bch)]))
        else:
            styoke[k] = scale_factor*sum(
                bch[0]['losses'][0]['stator']['stfe'][k])
            rotor[k] = scale_factor*sum(
                bch[0]['losses'][0]['rotor']['----'][k])

    if 'magnetH' in bch[0]['losses'][0]:
        k = 'magnetH'
    else:
        k = 'magnetJ'
    if len(pole_width) > 1:
        maglosses = _integrate1d(radius, scale_factor*np.array(
            [b['losses'][0][k]/l for l, b in zip(lfe, bch)]))
    else:
        maglosses = scale_factor*bch[0]['losses'][0][k]

    freq = bch[0]['losses'][0]['stator']['stfe']['freq'][0]

    wdg = windings.Winding(
        dict(Q=mmod.stator['num_slots'],
             p=mmod.poles//2,
             m=mmod.winding['num_phases'],
             l=mmod.winding['num_layers']))

    cufill = mmod.winding.get('cufilfact', 0.4)
    aw = (mmod.stator['afm_stator']['slot_width']*
              mmod.stator['afm_stator']['slot_height']*
              cufill/mmod.winding['num_wires']/mmod.winding['num_layers'])
    r1 = wdg_resistance(wdg, mmod.winding['num_wires'],
                        mmod.winding['num_par_wdgs'],
                        aw,
                        mmod.outer_diam, mmod.inner_diam)
    i1 = np.mean([np.max(c) for c in currents])/np.sqrt(2)
    plcu = mmod.winding['num_phases']*i1**2*r1

    return {
        'pos': pos.tolist(), 'r1': r1,
        'torque': torque,
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
        cu_losses = sum([b['losses'][0][wdgk] for b in bch])
        return scale_factor*cu_losses
    except KeyError:
        return 0  # noload calc has no winding losses


class AFPM:
    """Axial Flux PM
    Arguments:
    workdir: (str) working directory
    magnetizingCurves: (str) or (object)
    magnetMat: magnet material
    condMat: conductor material
    """
    def __init__(self, workdir, magnetizingCurves='.', magnetMat='',
                 condMat=''):
        self.parstudy = parstudy.List(
            workdir, condMat=condMat, magnets=magnetMat,
            magnetizingCurves=magnetizingCurves)

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
            for pw in pole_width:
                machine['num_agnodes'] = 6*round(pw/machine['airgap']/4)
                #Q = machine['stator']['num_slots']
                #p = machine['poles']
                #nper = np.lcm(Q, p)
                #if machine['num_agnodes'] < nper:
                #    machine['num_agnodes'] = 8*round(pw/machine['airgap']/4)

        parvardef = {
            "decision_vars": [
                {"values": pole_width,
                 "name": "pole_width"},
                {"values": lfe,
                 "name": "lfe"},
                {"values": linspeed, "name": "speed"}
            ]
        }
        machine['pole_width'] = np.pi * machine['inner_diam']/machine['poles']
        machine['lfe'] = machine['outer_diam'] - machine['inner_diam']

        nlresults = {}
        if (simulation['calculationMode'] != 'cogg_calc' and
            'poc' not in simulation):
            nlcalc = dict(
                calculationMode="cogg_calc",
                magn_temp=simulation.get('magn_temp', 20),
                num_move_steps=60,
                speed=0)
            logging.info("Noload simulation")
            nlresults = self.parstudy(parvardef,
                                      machine, nlcalc, engine)
            nlresults.update(process(lfe, pole_width, machine, nlresults['f']))

            current_angles = nlresults['f'][0]['current_angles']
            simulation['poc'] = poc.Poc(
                999,
                parameters={
                    'phi_voltage_winding': current_angles})
            logger.info("Current angles: %s", current_angles)

        lresults = self.parstudy(
            parvardef,
            {k: machine[k] for k in machine if k!= 'afm_rotor'},
            simulation, engine)  # Note: imcomplete machine prevents rebuild

        results = process(lfe, pole_width, machine, lresults['f'])
        results.update(_psidq_ldq(results, nlresults))
        return results
