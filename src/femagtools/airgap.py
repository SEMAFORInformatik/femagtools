# -*- coding: utf-8 -*-
"""
    femagtools.airgap
    ~~~~~~~~~~~~~~~~~

    Read airgap dat file


"""
import numpy as np
import logging
from . import utils

logger = logging.getLogger(__name__)


def fft(pos, b, pmod=0):
    """calculate fft spectrum of flux density and return samples,
    values, amplitude and phase of base harmonic

    Arguments:
      pos: (list of floats) sample positions
      b: (list of floats) flux density values
      pmod: number of poles in model (ignored if 0)
    """
    r = utils.fft(pos, b, pmod)
    Bamp = r['a']
    alfa0 = r['alfa0']
    T0 = r['T0']
    npoles = 2*int(np.ceil(360/T0))
    logger.info("flux density: %s poles B amp %f ",
                npoles, r['a'])
    return dict(Bamp=Bamp, npoles=npoles,
                phi0=alfa0,
                pos=pos.tolist(),
                B=b.tolist(),
                nue=np.arange(0, 9*npoles).tolist(),
                B_nue=r['nue'],
                B_fft=(Bamp*np.cos(2*np.pi*pos/T0+alfa0)).tolist(),
                Bi=r['yi'],
                phi=np.linspace(pos[0], 360+pos[0], len(r['yi'])).tolist())


def read(filename, pmod=0):
    """read dat file with columns (phi, Br, Bphi)
    returns samples, values, amplitude and phase of base harmonic

    Args:
      filename: the name of the file to be processed
      pmod: number of poles in model (ignored if 0)
    """
    bag = np.loadtxt(filename).T
    if len(bag) < 3:
        logger.warn("%s has incomplete content", filename)
        return(dict())

    return fft(bag[0], bag[1], pmod)
