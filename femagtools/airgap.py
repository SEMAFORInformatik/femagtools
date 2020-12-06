# -*- coding: utf-8 -*-
"""
    femagtools.airgap
    ~~~~~~~~~~~~~~~~~

    Read airgap dat file



"""
import numpy as np
import logging

logger = logging.getLogger(__name__)


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
    
    model_angle = bag[0][-1] - bag[0][0]
    ntiles = int(round(360/model_angle))

    if pmod:
        negative_periodic = pmod % 2
    else:
        negative_periodic = np.abs(np.sum(bag[1])/np.max(bag[1])) > 1

    if negative_periodic:
        bx = np.append(
            np.concatenate(
                [n*bag[1][:-1]
                 for n in [m % 2 or -1
                           for m in range(1, ntiles+1)]]),
            bag[1][0])
    else:
        bx = np.append(
            np.tile(bag[1][:-1], ntiles),
            bag[1][0])

    N = len(bx)
    # compute DFT from induction
    Y = np.fft.fft(bx)
    
    # find the peak (amplitude of base harmonic)
    i = np.argmax(np.abs(Y[:N//2]))
    a = 2*np.abs(Y[i])/N
    freq = np.fft.fftfreq(N, d=bag[0][1]-bag[0][0])
    T0 = np.abs(1/freq[i])
    npoles = 2*int(np.ceil(360/T0))
    logger.info("%s: %s poles B amp %f ",
                filename, npoles, a)

    alfa0 = np.angle(Y[i])
    return dict(Bamp=a, npoles=npoles,
                phi0=alfa0,
                pos=bag[0].tolist(),
                B=bag[1].tolist(),
                nue=np.arange(0, 9*npoles).tolist(),
                B_nue=(2*np.abs(Y[:9*npoles])/N).tolist(),
                B_fft=(a*np.cos(2*np.pi*bag[0]/T0+alfa0)).tolist(),
                Bi=bx.tolist(),
                phi=np.linspace(bag[0][0], 360+bag[0][0], len(bx)).tolist())
