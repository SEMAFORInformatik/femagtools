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
from .utils import betai1, iqd, invpark, K, T, puconv, dqpar_interpol
import copy
import logging

logger = logging.getLogger(__name__)


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


def create_from_eecpars(temp, eecpars, lfe=1, wdg=1):
    """create machine according to the eecpars:
    PM, EESM or IM"""
    rlfe = lfe
    rwdg = wdg
    if 'ldq' in eecpars:  # this is a PM (or EESM)
        if (isinstance(eecpars['ldq'], list) and
            'ex_current' in eecpars['ldq'][0] or
            isinstance(eecpars['ldq'], dict) and
                'ex_current' in eecpars['ldq']):
            pars = copy.deepcopy(eecpars)
            pars['tcu1'] = temp[0]
            pars['tcu2'] = temp[1]
            machine = femagtools.machine.sm.SynchronousMachine(pars, lfe=rlfe, wdg=rwdg)
            try:
                machine.rotor_mass = rlfe*eecpars['rotor_mass']
            except KeyError:
                pass
            return machine

        if isinstance(eecpars['ldq'], list) and len(eecpars['ldq']) > 1:
            x, dqp = dqpar_interpol(
                temp[1], eecpars['ldq'], ipkey='temperature')
        else:
            dqp = eecpars['ldq'][0]
            logger.warning(
                "single temperature DQ parameters: unable to fit temperature %s", temp)

        beta = dqp['beta']
        i1 = np.array(dqp['i1'])/rwdg
        psid = rwdg*rlfe*dqp['psid']
        psiq = rwdg*rlfe*dqp['psiq']
        try:
            losses = __scale_losses(dqp['losses'], rlfe)
            losses['ef'] = eecpars['ldq'][-1]['losses']['ef']
            losses['eh'] = eecpars['ldq'][-1]['losses']['ef']
        except KeyError as e:
            logger.warning(e)
            losses = {}

        machine = PmRelMachineLdq(
            eecpars['m'], eecpars['p'],
            r1=eecpars['r1']*rlfe*rwdg**2,
            ls=eecpars['ls1']*rwdg**2,
            psid=psid,
            psiq=psiq,
            losses=losses,
            beta=beta,
            i1=i1,
            tcu1=temp[0])
        try:
            machine.rotor_mass = rlfe*eecpars['rotor_mass']
        except KeyError:
            pass
        return machine

    # must be an induction machine (TODO: check scaling)
    pars = copy.deepcopy(eecpars)
    pars['r1'] = rlfe*rwdg**2*pars['r1']
    pars['lsigma1'] = rlfe*pars['lsigma1']
    pars['lsigma2'] = rlfe*pars['lsigma2']
    pars['psiref'] = rwdg*rlfe*pars['psiref']
    pars['u1ref'] = rwdg*rlfe*pars['u1ref']
    pars['rotor_mass'] = rlfe*pars['rotor_mass']
    pars['r2'] = rlfe*pars['r2']
    pars['fec'] = rlfe*pars['fec']
    pars['fee'] = rlfe*pars['fee']
    pars['im'] = [im/rwdg for im in pars['im']]
    pars['psi'] = [psi*rwdg*rlfe for psi in pars['psi']]
    pars['tcu1'] = temp[0]
    pars['tcu2'] = temp[1]
    return InductionMachine(pars)


def __scale_losses(losses, rlfe):
    if losses:
        l = {k: rlfe*np.array(losses[k]) for k in (
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
      rlfe: scale factor length
      rwdg: scale factor number of windings
"""
    rwdg = wdg
    rlfe = lfe
    m = 3
    if isinstance(bch, Reader):
        p = bch.machine['p']
        if bch.type.lower().find('psid-psiq-identification') >= 0:
            id = np.array(bch.psidq['id'])/rwdg
            iq = np.array(bch.psidq['iq'])/rwdg
            psid = rwdg*rlfe*np.array(bch.psidq['psid'])
            psiq = rwdg*rlfe*np.array(bch.psidq['psiq'])
            try:
                losses = __scale_losses(bch.psidq['losses'], rlfe)
                losses['ef'] = bch.lossPar.get('ef', [2.0, 2.0])
                losses['eh'] = bch.lossPar.get('eh', [1.0, 1.0])
            except KeyError:
                losses = {}
            machine = PmRelMachinePsidq(m, p, psid, psiq, r1*rlfe*rwdg**2,
                                        id, iq, ls*rwdg**2, losses=losses)
            try:
                machine.rotor_mass = rlfe*np.sum(bch.weights[1])
            except IndexError:
                pass
            return machine


        if bch.type.lower().find('ld-lq-identification') >= 0:
            beta = bch.ldq['beta']
            i1 = np.array(bch.ldq['i1'])/rwdg
            psid = rwdg*rlfe*np.array(bch.ldq['psid'])
            psiq = rwdg*rlfe*np.array(bch.ldq['psiq'])
            try:
                losses = __scale_losses(bch.ldq['losses'], rlfe)
                losses['ef'] = bch.lossPar.get('ef', [2.0, 2.0])
                losses['eh'] = bch.lossPar.get('eh', [1.0, 1.0])
            except KeyError:
                losses = {}
            machine = PmRelMachineLdq(m, p, psid=psid, psiq=psiq,
                                      r1=r1*rlfe*rwdg**2,
                                      i1=i1, beta=beta, ls=ls*rwdg**2,
                                      losses=losses)
            try:
                machine.rotor_mass = rlfe*np.sum(bch.weights[1])
            except IndexError:
                pass
            return machine

        raise ValueError("Unsupported BCH type {}".format(bch.type))
    # must be ERG type:
    p = int(round(np.sqrt(2)*bch['M_sim'][-1][-1]/(
        m*bch['Psi_d'][-1][-1] * bch['i1'][-1])))

    machine = PmRelMachineLdq(m, p, r1=r1*rlfe*rwdg**2,
                              beta=bch['beta'], i1=np.array(bch['i1'])/rwdg,
                              psid=rwdg*rlfe*np.array(bch['Psi_d'])/np.sqrt(2),
                              psiq=rwdg*rlfe*np.array(bch['Psi_q'])/np.sqrt(2),
                              ls=ls*rwdg**2)
    try:
        machine.rotor_mass = rlfe*np.sum(bch.weights[1])
    except IndexError:
        pass
    return machine
