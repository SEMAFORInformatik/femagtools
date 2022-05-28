# -*- coding: utf-8 -*-
"""
    femagtools.machine
    ~~~~~~~~~~~~~~~~~~

    Analytical machine models

"""
import numpy as np
from femagtools.bch import Reader
from .pm import PmRelMachineLdq, PmRelMachinePsidq
from .utils import betai1, iqd, invpark, K, T, puconv, dqpar_interpol, wdg_resistance


def __scale_losses(losses, wdg, lfe):
    if losses:
        l = {k: wdg*lfe*np.array(losses[k]) for k in (
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
                losses = __scale_losses(bch.psidq['losses'], wdg, lfe)
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
                losses = __scale_losses(bch.ldq['losses'], wdg, lfe)
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
