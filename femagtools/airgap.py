# -*- coding: utf-8 -*-
"""
    femagtools.airgap
    ~~~~~~~~~~~~~~~~~

    Read airgap dat file



"""
import numpy as np
import logging

logger = logging.getLogger(__name__)


def read(filename, pmod):
    """read dat file with columns (phi, Br, Bphi)
    returns samples, values, amplitude and phase of base harmonic

    Args:
      filename: the name of the file to be processed
      pmod: number of poles in model
    """
    bag = np.loadtxt(filename).T
    if len(bag) < 3:
        logger.warn("%s has incomplete content", filename)
        return(dict())

    N = 2**10  # The DFT is most efficient when N is a power of 2
    phi = np.linspace(0, 2*np.pi, N)
    nphi = int(360/(bag[0][1] - bag[0][0])) + 1
    npoles = (nphi-1)//(len(bag[0]) - 1)
    phitab = np.linspace(0, 2*np.pi, nphi)
    if pmod % 2:
        bx = np.append(
            np.concatenate(
                [n*bag[1][:-1]
                 for n in [m % 2 or -1
                           for m in range(1, npoles+1)]]),
            bag[1][0])
    else:
        bx = np.append(
            np.tile(bag[1][:-1], npoles),
            bag[1][0])

    br = np.interp(phi, phitab, bx)
        
    # compute DFT from induction
    Y = np.fft.fft(br)

    # find the peak (amplitude of base harmonic)
    i = np.argmax(np.abs(Y))
    a = 2*np.abs(Y[i])/N
    logger.info("%s: B amp %f",
                filename, a)
    
    alfa0 = np.angle(Y[i])
    alfa = npoles/2*bag[0]/180*np.pi

    return dict(Bamp=a,
                phi0=alfa0,
                pos=bag[0],
                B=bag[1],
                nue=np.arange(0, 9*npoles),
                B_nue=2*np.abs(Y[:9*npoles])/N,
                B_fft=(a*np.cos(alfa+alfa0)).tolist())
