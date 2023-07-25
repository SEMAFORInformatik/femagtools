"""
    femagtools.machine.afpm
    ~~~~~~~~~~~~~~~~~~~~~~~

    Axial Flux PM Machine

"""
import logging
import numpy as np
from .. import poc
from .. import parstudy
from .. import model
from .. import utils
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

def integrate(radius, pos, val):
    interp = RegularGridInterpolator((radius, pos), val)
    def func(x, y):
        return interp((x, y))
    return [quad(func, radius[0], radius[-1], args=(p,))[0]
            for p in pos]

def integrate1d(radius, val):
    interp = interp1d(radius, val)
    def func(x):
        return interp((x))
    return quad(func, radius[0], radius[-1])[0]


def process(lfe, pole_width, machine, bch):
    # process results: torque, voltage, losses
    model_type = machine['afmtype']
    num_slots = machine['stator']['num_slots']
    mmod = model.MachineModel(machine)
    slots_gen = mmod.stator['num_slots_gen']
    scale_factor = get_scale_factor(model_type, num_slots, slots_gen)
    endpos = [2*pw*1e3 for pw in pole_width]
    displ = [[d for d in r['linearForce'][0]['displ']
              if d < e*(1+1/len(r['linearForce'][0]['displ']))]
             for e, r in zip(endpos, bch)]
    radius = [pw*machine['poles']/2/np.pi for pw in pole_width]
    rotpos = [np.array(d)/r/1000 for r, d in zip(radius, displ)]
    torque = [r*scale_factor*np.array(fx)/l
              for l, r, fx in zip(lfe, radius,
                                  [r['linearForce'][0]['force_x']
                                   for r in bch])]
    voltage = {k: [scale_factor * np.array(ux)/l
                   for l, ux in zip(lfe, [r['flux'][k][0]['voltage_dpsi']
                                          for r in bch])]
               for k in bch[0]['flux']}

    n = len(rotpos[0])
    currents = [bch[0]['flux'][k][0]['current_k'][:n]
                for k in bch[0]['flux']]
    emf = [integrate(radius, rotpos[0], np.array(voltage[k])[:, :n])
           for k in voltage]
    pos = (rotpos[0]/np.pi*180)
    emffft = utils.fft(pos, np.array(emf[0]))

    pfelosses = [integrate1d(radius, scale_factor*np.array(
        [(b['losses'][0]["staza"] + b['losses'][0]["stajo"])/l
         for l, b in zip(lfe, bch)])),
                 integrate1d(radius, scale_factor*np.array(
                     [b['losses'][0]['rotfe']/l
                      for l, b in zip(lfe, bch)]))]
    maglosses = integrate1d(radius, scale_factor*np.array(
        [b['losses'][0]["magnetJ"]/l for l, b in zip(lfe, bch)]))

    freq = bch[0]['losses'][0]['stator']['stfe']['freq'][0]

    return {
        'pos': pos.tolist(),
        'torque': integrate(radius, rotpos[0], np.array(torque)[:, :n]),
        'emf': emf,
        'emf_amp': emffft['a'], 'emf_angle': emffft['alfa0'],
        'freq': freq,
        'currents': currents,
        'plfe': pfelosses, 'plmag': maglosses,
        'plcu': get_copper_losses(scale_factor, bch)}


def get_scale_factor(model_type, num_slots, slots_gen):
    """Determines the scale factor
    Parameters
    ----------
    model_type : (str) type of model
    num_slots: (int) total number of stator slots
    slots_gen: (int) number of slots in model

    Return
    ------
    scale_factor : int
        Scale factor based on number of poles, number of poles simulated
        and the model type
    """
    segments = num_slots / slots_gen
    if model_type in ("S2R1", "S1R2"):
        return 2 * segments
    return segments


def get_arm_lengths(outer_diam, inner_diam, num_slices):
    d = outer_diam - inner_diam
    return [d/(4*(num_slices-1))] + [
        d/(2*(num_slices-1))
        for i in range(1,num_slices-1)] + [d/(4*(num_slices-1))]


def get_pole_widths(outer_diam, inner_diam, poles, num_slices):
    d = outer_diam - inner_diam
    return [np.pi * inner_diam/poles] + [
        np.pi * (inner_diam + d*i/(num_slices - 1))/poles
        for i in range(1, num_slices-1)] + [np.pi * outer_diam/poles]


def get_copper_losses(scale_factor, bch):
    # relative length of winding = 1
    # winding losses are without winding head
    try:
        cu_losses = sum([b['losses'][0]['winding'] for b in bch])
        return scale_factor*cu_losses
    except KeyError:
        return 0  # noload calc has no winding losses


class AFPM:
    def __init__(self, workdir, magnetizingCurves='.', magnetMat='',
                 condMat=''):
        self.parstudy = parstudy.List(
            workdir, condMat=condMat, magnets=magnetMat,
            magnetizingCurves=magnetizingCurves)

    def __call__(self, engine, machine, simulation, num_slices):
        # check afm type
        if machine['afmtype'] not in AFM_TYPES:
            raise ValueError(f"invalid afm type {machine['afmtype']}")

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
                machine['num_agnodes'] = 3*round(pw/machine['airgap']/2)

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

        if simulation['calculationMode'] != 'cogg_calc':
            nlcalc = dict(
                calculationMode="cogg_calc",
                magn_temp=simulation.get('magn_temp', 20),
                num_move_steps=60,
                speed=0)
            logging.info("Noload simulation")
            nlresults = self.parstudy(parvardef, machine, nlcalc, engine)
            nlresults.update(process(lfe, pole_width, machine, nlresults['f']))

            if 'poc' not in simulation:
                current_angles = nlresults['f'][0]['current_angles']
                simulation['poc'] = poc.Poc(
                    999,
                    parameters={
                        'phi_voltage_winding': current_angles})
                logger.info("Current angles: %s", current_angles)

        results = self.parstudy(parvardef, machine, simulation, engine)
        results.update(process(lfe, pole_width, machine, results['f']))
        gamma = -(results['emf_angle'] - nlresults['emf_angle'])
        w1 = 2*np.pi*results['freq']
        results['psid'] = np.cos(gamma)*results['emf_amp']/w1
        results['psiq'] = np.sin(gamma)*results['emf_amp']/w1
        results['psim'] = nlresults['emf_amp']/w1
        results['i1'] = np.mean([np.max(c)
                                 for c in results['currents']])/np.sqrt(2)
        beta = results['f'][0]['losses'][0]['beta']/180*np.pi
        results['beta'] = beta/np.pi*180
        results['id'] = np.sqrt(2)*results['i1']*np.cos(beta)
        results['iq'] = np.sqrt(2)*results['i1']*np.sin(beta)
        if np.abs(results['id']) > 0:
            results['Ld'] = (results['psid'] - results['psim'])/results['id']
        if np.abs(results['iq']) > 0:
            results['Lq'] = results['psiq']/results['iq']
        return results
