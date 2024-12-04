""":mod:`femagtools.machine` -- Analytical machine models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
import numpy as np
from femagtools.bch import Reader
from .pm import PmRelMachineLdq, PmRelMachinePsidq
from .sm import SynchronousMachine, SynchronousMachineLdq, SynchronousMachinePsidq
from .im import InductionMachine
from .utils import betai1, iqd, invpark, K, T, puconv, dqpar_interpol
import copy
import logging

logger = logging.getLogger(__name__)


def __scale_losses(losses, rlfe):
    if losses:
        l = {k: rlfe*np.array(losses[k]) for k in (
            'styoke_hyst', 'styoke_eddy',
            'stteeth_hyst', 'stteeth_eddy',
            'styoke_excess', 'stteeth_excess', 'rotor_excess',
            'rotor_hyst', 'rotor_eddy',
            'magnet') if k in losses}
        l['speed'] = losses['speed']
        return l
    return {}


def create_from_eecpars(temp, eecpars, lfe=1, wdg=1):
    """create machine according to the eecpars:
    PM, EESM or IM"""
    rlfe = lfe
    rwdg = wdg
    opts = {k: eecpars[k] for k in ('zeta1', 'gam', 'kh', 'kpfe',
                                    'kfric_b', 'kpmag') if k in eecpars}
    try:
        opts['rotor_mass'] = rlfe*eecpars['rotor_mass']
    except KeyError:
        pass

    if 'ldq' in eecpars or 'psidq' in eecpars:  # this is a PM (or EESM)
        try:
            dqpars = eecpars['ldq']
        except KeyError:
            dqpars = eecpars['psidq']
        if (isinstance(dqpars, list) and
            'ex_current' in dqpars[0] or
            isinstance(dqpars, dict) and
                'ex_current' in dqpars):
            smpars = copy.deepcopy(eecpars)
            smpars['tcu1'] = temp[0]
            smpars['tcu2'] = temp[1]
            if 'ldq' in smpars:
                machine = SynchronousMachineLdq(smpars, lfe=rlfe, wdg=rwdg, **opts)
            else:
                machine = SynchronousMachinePsidq(smpars, lfe=rlfe, wdg=rwdg, **opts)
            return machine

        if isinstance(dqpars, list) and len(dqpars) > 1:
            x, dqp = dqpar_interpol(
                temp[1], dqpars, ipkey='temperature')
        else:
            dqp = dqpars[0]
            logger.warning(
                "single temperature DQ parameters: unable to fit temperature %s", temp)

        psid = rwdg*rlfe*dqp['psid']
        psiq = rwdg*rlfe*dqp['psiq']
        losses = __scale_losses(dqp['losses'], rlfe)
        losses['ef'] = dqpars[-1]['losses'].get('ef', [2.0, 2.0])
        losses['hf'] = dqpars[-1]['losses'].get('hf', [1.0, 1.0])
        # TODO handle bertotti excess loss factor

        if 'psidq' in eecpars:
            machine = PmRelMachinePsidq(
                eecpars['m'], eecpars['p'],
                r1=eecpars.get('r1', 0)*rlfe*rwdg**2,
                ls=eecpars.get('ls1', 0)*rwdg**2,
                psid=psid,
                psiq=psiq,
                losses=losses,
                id=np.array(dqp['id'])/rwdg,
                iq=np.array(dqp['iq'])/rwdg,
                tcu1=temp[0],
                **opts)
        else:
            beta = dqp['beta']
            i1 = np.array(dqp['i1'])/rwdg
            machine = PmRelMachineLdq(
                eecpars['m'], eecpars['p'],
                r1=eecpars.get('r1', 0)*rlfe*rwdg**2,
                ls=eecpars.get('ls1', 0)*rwdg**2,
                psid=psid,
                psiq=psiq,
                losses=losses,
                beta=beta,
                i1=i1,
                tcu1=temp[0],
                **opts)
        return machine

    # must be an induction machine (TODO: check scaling)
    pars = copy.deepcopy(eecpars)
    pars['r1'] = rlfe*rwdg**2*pars.get('r1', 0)
    pars['lsigma1'] = rlfe*pars['lsigma1']
    pars['lsigma2'] = rlfe*pars['lsigma2']
    pars['psiref'] = rwdg*rlfe*pars['psiref']
    pars['u1ref'] = rwdg*rlfe*pars['u1ref']
    pars['r2'] = rlfe*pars['r2']
    pars['fec'] = rlfe*pars['fec']
    pars['fee'] = rlfe*pars['fee']
    pars['im'] = [im/rwdg for im in pars['im']]
    pars['psi'] = [psi*rwdg*rlfe for psi in pars['psi']]
    pars['tcu1'] = temp[0]
    pars['tcu2'] = temp[1]
    pars.update(opts)
    return InductionMachine(pars)


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
                losses['hf'] = bch.lossPar.get('hf', [1.0, 1.0])
            except KeyError:
                losses = {}
            if 'ex_current' in bch.machine:
                raise ValueError("not yet implemented for EESM")
            machine = PmRelMachinePsidq(m, p, psid, psiq, r1*rlfe*rwdg**2,
                                        id, iq, ls*rwdg**2, losses=losses)
            try:
                machine.rotor_mass = rlfe*np.sum(bch.weights[1])
            except (IndexError, AttributeError):
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
                losses['hf'] = bch.lossPar.get('hf', [1.0, 1.0])
            except KeyError:
                losses = {}
            if 'ex_current' in bch.machine:
                raise ValueError("not yet implemented for EESM")

            machine = PmRelMachineLdq(m, p, psid=psid, psiq=psiq,
                                      r1=r1*rlfe*rwdg**2,
                                      i1=i1, beta=beta, ls=ls*rwdg**2,
                                      losses=losses)
            try:
                machine.rotor_mass = rlfe*np.sum(bch.weights[1])
            except (IndexError, AttributeError):
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
    except (IndexError, AttributeError):
        pass
    return machine
