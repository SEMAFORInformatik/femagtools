# -*- coding: utf-8 -*-
"""
    femagtools.machine
    ~~~~~~~~~~~~~~~~~~

    Analytical machine models

"""
import numpy as np
from femagtools.bch import Reader
from .pm import PmRelMachineLdq, PmRelMachinePsidq
from .im import InductionMachine
from .utils import betai1, iqd, invpark, K, T, puconv, dqpar_interpol, wdg_resistance
import femagtools.parstudy
import copy


def create_from_eecpars(temp, eecpars):
    """create machine according to the eecpars:
    PM, EESM or IM"""
    if 'ldq' in eecpars:  # this is a PM (or EESM)
        if (isinstance(eecpars['ldq'], list) and
            'ex_current' in eecpars['ldq'][0] or
            isinstance(eecpars['ldq'], dict) and
                'ex_current' in eecpars['ldq']):
            pars = copy.deepcopy(eecpars)
            pars['tcu1'] = temp[0]
            pars['tcu2'] = temp[1]
            return femagtools.machine.sm.SynchronousMachine(pars)

        if isinstance(eecpars['ldq'], list) and len(eecpars['ldq']) > 1:
            x, dqp = dqpar_interpol(
                temp[1], eecpars['ldq'], ipkey='temperature')
        else:
            dqp = eecpars['ldq'][0]
            logger.warning(
                "single temperature DQ parameters: unable to fit temperature %s", temp)
        return PmRelMachineLdq(
            eecpars['m'], eecpars['p'], r1=eecpars['r1'],
            ls=eecpars['ls1'],
            psid=dqp['psid'],
            psiq=dqp['psiq'],
            losses=dqp['losses'],
            beta=dqp['beta'],
            i1=dqp['i1'],
            tcu1=temp[0])
    else:  # must be an induction machine
        pars = copy.deepcopy(eecpars)
        pars['tcu1'] = temp[0]
        pars['tcu2'] = temp[1]
        return InductionMachine(pars)


def __scale_losses(losses, lfe):
    if losses:
        l = {k: lfe*np.array(losses[k]) for k in (
            'styoke_hyst', 'styoke_eddy',
            'stteeth_hyst', 'stteeth_eddy',
            'rotor_hyst', 'rotor_eddy',
            'magnet')}
        l['speed'] = losses['speed']
        return l
    return {}


def create(bch, r1, ls, lfe=1, wdg=1):
    """create PmRelMachine from BCH

    Arguments:
      bch: BchReader or Erg object
      r1: winding resistance
      ls: winding leakage
      lfe: scale factor length
      wdg: scale factor number of windings
"""
    m = 3
    if isinstance(bch, Reader):
        p = bch.machine['p']
        if bch.type.lower().find('psid-psiq-identification') >= 0:
            id = np.array(bch.psidq['id'])/wdg
            iq = np.array(bch.psidq['iq'])/wdg
            psid = wdg*lfe*np.array(bch.psidq['psid'])
            psiq = wdg*lfe*np.array(bch.psidq['psiq'])
            try:
                losses = __scale_losses(bch.psidq['losses'], lfe)
                losses['ef'] = bch.lossPar.get('ef', [2.0, 2.0])
                losses['eh'] = bch.lossPar.get('eh', [1.0, 1.0])
            except KeyError:
                losses = {}
            return PmRelMachinePsidq(m, p, psid, psiq, r1*lfe*wdg**2,
                                     id, iq, ls*wdg**2, losses=losses)

        if bch.type.lower().find('ld-lq-identification') >= 0:
            beta = bch.ldq['beta']
            i1 = np.array(bch.ldq['i1'])/wdg
            psid = wdg*lfe*np.array(bch.ldq['psid'])
            psiq = wdg*lfe*np.array(bch.ldq['psiq'])
            try:
                losses = __scale_losses(bch.ldq['losses'], lfe)
                losses['ef'] = bch.lossPar.get('ef', [2.0, 2.0])
                losses['eh'] = bch.lossPar.get('eh', [1.0, 1.0])
            except KeyError:
                losses = {}
            return PmRelMachineLdq(m, p, psid=psid, psiq=psiq,
                                   r1=r1*lfe*wdg**2,
                                   i1=i1, beta=beta, ls=ls*wdg**2,
                                   losses=losses)
        raise ValueError("Unsupported BCH type {}".format(bch.type))
    # must be ERG type:
    p = int(round(np.sqrt(2)*bch['M_sim'][-1][-1]/(
        m*bch['Psi_d'][-1][-1] * bch['i1'][-1])))

    return PmRelMachineLdq(m, p, r1=r1*lfe*wdg**2,
                           beta=bch['beta'], i1=np.array(bch['i1'])/wdg,
                           psid=wdg*lfe*np.array(bch['Psi_d'])/np.sqrt(2),
                           psiq=wdg*lfe*np.array(bch['Psi_q'])/np.sqrt(2),
                           ls=ls*wdg**2)


def dqparident(workdir, engine, temp, machine,
               magnetizingCurves, magnetMat=[], condMat=[],
               **kwargs):
    """return list of parameters of equivalent circuit for PM machines

    arguments:
    workdir -- directory for intermediate files
    engine -- calculation driver (multiproc, docker, condor)

    temp -- list of magnet temperatures in degree Celsius
    machine -- dict() with machine parameters
    magnetizingCurves -- list of dict() with BH curves
    magnetMat -- list of dict() with magnet material properties
    condMat -- list of dict() with conductor material properties

    optional arguments:
    num_cur_steps: number of current steps (default 5)
    num_beta_steps: number of current steps (default 13)
    speed: rotor speed in 1/s (default 160/p)
    i1_max: maximum current in A rms (default approx 3*i1nom)
    use_multiprocessing: (boolean) perfom FE simulations in parallel (default: True)
    """
    import pathlib

    da1 = machine['outer_diam']
    Q1 = machine['stator']['num_slots']
    if 'statorRotor3' in machine['stator']:
        hs = machine['stator']['statorRotor3']['slot_height']
    elif 'stator1' in machine['stator']:
        hs = machine['stator']['stator1']['slot_rf1'] - \
            machine['stator']['stator1']['tip_rh1']
    elif 'stator4' in machine['stator']:
        hs = machine['stator']['stator4']['slot_height']
    N = machine['windings']['num_wires']
    Jmax = 15  # max current density in A/mm2

    i1_max = round(0.28*np.pi*hs*(da1+hs)/Q1/N*Jmax*1e5)*10 * \
        machine['windings'].get('num_par_wdgs', 1)
    period_frac = 6
    if machine.get('external_rotor', False):
        period_frac = 1  # TODO: missing femag support

    # winding resistance
    yd = machine['windings'].get('coil_span', Q1/machine['poles'])
    wdg = femagtools.windings.Winding(
        {'Q': machine['stator']['num_slots'],
         'm': machine['windings']['num_phases'],
         'p': machine['poles']//2,
         'l': machine['windings']['num_layers'],
         'yd': yd})

    lfe = machine['lfe']
    g = machine['windings'].get('num_par_wdgs', 1)
    aw = np.pi*machine['windings'].get('dia_wire', 1e-3)**2/4
    r1 = wdg_resistance(wdg, N, g, aw, da1, hs, lfe)

    parvardef = {
        "decision_vars": [
            {"values": temp, "name": "magn_temp"},
        ]
    }

    parvar = femagtools.parstudy.List(
        workdir, condMat=condMat,
        magnetizingCurves=magnetizingCurves,
        magnets=magnetMat)

    simulation = dict(
        calculationMode='ld_lq_fast',
        i1_max=kwargs.get('i1_max', i1_max),
        magn_temp=20,
        wind_temp=20,
        beta_max=0.0,
        beta_min=-180.0,
        num_par_wdgs=machine['windings'].get('num_par_wdgs', 1),
        num_cur_steps=kwargs.get('num_cur_steps', 5),
        num_beta_steps=kwargs.get('num_beta_steps', 13),
        speed=kwargs.get('speed', 160/machine['poles']),
        period_frac=period_frac)

    # TODO: cleanup()  # remove previously created files in workdir
    # start calculation
    results = parvar(parvardef, machine, simulation, engine)
    import json
    with open('results.json', 'w') as fp:
        json.dump(results, fp)
    ls1 = 0
    leakfile = pathlib.Path(workdir) / 'end_wind_leak.dat'
    if leakfile.exists():
        leakages = leakfile.read_text().split()
        ls1 = np.linalg.norm(leakages[1:])
    try:
        rotor_mass = sum(results['f'][-1]['weights'][-1])
    except KeyError:
        rotor_mass = 0  # need femag classic > rel-9.3.x-48-gca42bbd0

    ldq = [dict(
        i1=b['ldq']['i1'],
        beta=b['ldq']['beta'],
        psid=b['ldq']['psid'],
        psiq=b['ldq']['psiq'],
        torque=b['ldq']['torque'],
        ld=b['ldq']['ld'],
        lq=b['ldq']['lq'],
        psim=b['ldq']['psim']) for b in results['f']]

    losskeys = ('speed',
                'styoke_hyst', 'stteeth_hyst', 'styoke_eddy',
                'stteeth_eddy', 'rotor_hyst', 'rotor_eddy',
                'magnet')
    for i in range(len(results['f'])):
        ldq[i]['temperature'] = results['x'][0][i]
        ldq[i]['losses'] = {k: results['f'][i]['ldq']['losses'][k]
                            for k in losskeys}
        for k in ('hf', 'ef'):
            ldq[i]['losses'][k] = results['f'][i]['lossPar'][k]

    return {'m': machine['windings']['num_phases'],
            'p': machine['poles']//2,
            'r1': machine['windings'].get('resistance', r1),
            'ls1': ls1,
            "rotor_mass": rotor_mass, "kfric_b": 1,
            'ldq': ldq}
