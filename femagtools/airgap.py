# -*- coding: utf-8 -*-
"""
    femagtools.airgap
    ~~~~~~~~~~~~~~~~~

    Read airgap dat file


"""
import numpy as np
import logging

logger = logging.getLogger(__name__)


def fft(pos, b, pmod=0):
    """calculate fft spectrum of flux density and return samples,
    values, amplitude and phase of base harmonic

    Arguments:
      pos: (list of floats) sample positions
      b: (list of floats) flux density values
      pmod: number of poles in model (ignored if 0)
    """
    model_angle = pos[-1] - pos[0]
    ntiles = int(round(360/model_angle))

    if pmod:
        negative_periodic = pmod % 2
    else:
        negative_periodic = np.abs(b[0] - b[-1])/np.max(b) > 1

    if negative_periodic:
        bx = np.append(
            np.concatenate(
                [n*b[:-1]
                 for n in [m % 2 or -1
                           for m in range(1, ntiles+1)]]),
            b[0])
    else:
        bx = np.append(
            np.tile(b[:-1], ntiles),
            b[0])

    N = len(bx)
    # compute DFT from induction
    Y = np.fft.fft(bx)

    # find the peak (amplitude of base harmonic)
    i = np.argmax(np.abs(Y[:N//2]))
    a = 2*np.abs(Y[i])/N
    freq = np.fft.fftfreq(N, d=pos[1]-pos[0])
    T0 = np.abs(1/freq[i])
    npoles = 2*int(np.ceil(360/T0))
    logger.info("flux density: %s poles B amp %f ",
                npoles, a)

    alfa0 = np.angle(Y[i])
    return dict(Bamp=a, npoles=npoles,
                phi0=alfa0,
                pos=pos.tolist(),
                B=b.tolist(),
                nue=np.arange(0, 9*npoles).tolist(),
                B_nue=(2*np.abs(Y[:9*npoles])/N).tolist(),
                B_fft=(a*np.cos(2*np.pi*pos/T0+alfa0)).tolist(),
                Bi=bx.tolist(),
                phi=np.linspace(pos[0], 360+pos[0], len(bx)).tolist())


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
