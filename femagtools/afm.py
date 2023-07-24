"""
    femagtools.afm
    ~~~~~~~~~~~~~~~~

    Axial Flux Machine

"""
import logging
import numpy as np
import femagtools.parstudy
from scipy.interpolate import RegularGridInterpolator
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


def process(lfe, pole_width, machine, results):
    # process results: torque, voltage, losses
    model_type = machine['afmtype']
    num_slots = machine['stator']['num_slots']
    model = femagtools.model.MachineModel(machine)
    slots_gen = model.stator['num_slots_gen']
    scale_factor = get_scale_factor(model_type, num_slots, slots_gen)
    endpos = [2*pw*1e3 for pw in pole_width]
    displ = [[d for d in r['linearForce'][0]['displ']
              if d < e*(1+1/len(r['linearForce'][0]['displ']))]
             for e, r in zip(endpos, results['f'])]
    radius = [pw*machine['poles']/2/np.pi for pw in pole_width]
    rotpos = [np.array(d)/r/1000 for r, d in zip(radius, displ)]
    torque = [r*scale_factor*np.array(fx)/l
              for l, r, fx in zip(lfe, radius,
                                  [r['linearForce'][0]['force_x']
                                   for r in results['f']])]
    voltage = {k: [scale_factor * np.array(ux)/l
                   for l, ux in zip(lfe, [r['flux'][k][0]['voltage_dpsi']
                                          for r in results['f']])]
               for k in results['f'][0]['flux']}
    n = len(rotpos[0])
    return {
        'pos': (rotpos[0]/np.pi*180).tolist(),
        'torque': integrate(radius, rotpos[0], np.array(torque)[:, :n]),
        'voltage': [integrate(radius, rotpos[0], np.array(voltage[k])[:, :n])
                    for k in voltage]}


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


def get_arm_lengths(outer_diam, inner_diam, poles, num_slices):
    d = outer_diam - inner_diam
    return [d/(4*(num_slices-1))] + [
        d/(2*(num_slices-1))
        for i in range(1,num_slices-1)] + [d/(4*(num_slices-1))]


def get_pole_widths(outer_diam, inner_diam, poles, num_slices):
    d = outer_diam - inner_diam
    return [np.pi * inner_diam/poles] + [
        np.pi * (inner_diam + d*i/(num_slices - 1))/poles
        for i in range(1, num_slices-1)] + [np.pi * outer_diam/poles]


class AFM:
    def __init__(self, workdir, magnetizingCurves='.', magnetMat='',
                 condMat=''):
        self.parstudy = femagtools.parstudy.List(
            workdir, condMat=condMat, magnets=magnetMat,
            magnetizingCurves=magnetizingCurves)

    def cogg_calc(self, speed, machine):
        machine['pole_width'] = np.pi * machine['inner_diam']/machine['poles']
        machine['lfe'] = machine['outer_diam'] - machine['inner_diam']

        simulation = dict(
            calculationMode="cogg_calc",
            magn_temp=20.0,
            num_move_steps=60,
            speed=speed*machine['pole_width']*machine['poles'])
        logging.info("Noload simulation")
        return self.parstudy.femag(machine, simulation)

    def __call__(self, engine, machine, simulation, num_slices):
        # check afm type
        if machine['afmtype'] not in AFM_TYPES:
            raise ValueError(f"invalid afm type {machine['afmtype']}")

        if 'poc' not in simulation:
            rcogg = self.cogg_calc(simulation['speed'], machine)
            simulation['poc'] = femagtools.poc.Poc(
                999,
                parameters={
                    'phi_voltage_winding': rcogg.current_angles})
            logger.info("Current angles: %s", rcogg.current_angles)

        lfe = get_arm_lengths(machine['outer_diam'],
                              machine['inner_diam'],
                              machine['poles'],
                              num_slices)
        pole_width = get_pole_widths(machine['outer_diam'],
                                     machine['inner_diam'],
                                     machine['poles'],
                                     num_slices)
        linspeed = [simulation['speed']*machine['poles']*pw
                    for pw in pole_width]

        parvardef = {
            "decision_vars": [
                {"values": pole_width,
                 "name": "pole_width"},
                {"values": lfe,
                 "name": "lfe"},
                {"values": linspeed, "name": "speed"}
            ]
        }
        results = self.parstudy(parvardef, machine, simulation, engine)
        results.update(process(lfe, pole_width, machine, results))
        return results
