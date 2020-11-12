# -*- coding: utf-8 -*-
"""
    femagtools.airgap
    ~~~~~~~~~~~~~~~~~

    Read airgap dat file



"""
import numpy as np
import logging
import math
import time

logger = logging.getLogger(__name__)


def read(filename, pmod):
    """read dat file with columns (phi, Br, Bphi)
    returns samples, values, amplitude and phase of base harmonic

    Args:
      filename: the name of the file to be processed
      pmod:     number of poles in model
    """
    bag = np.loadtxt(filename).T
    if len(bag) < 3:
        logger.warn("%s has incomplete content", filename)
        return(dict())

    time1 = time.time()
    model_angle = bag[0][-1]-bag[0][0]
    pole_pitch = model_angle/pmod
    pole_pairs = int(math.ceil(360.0/pole_pitch-0.5)/2)
    per_bnd = math.ceil(model_angle/pole_pitch-0.5) % 2
    n_copy = math.ceil(360.0/model_angle-0.5)
    
    bi = []    
    xi = []
    if per_bnd:
        ''' 
          Periodicity of the model:
          (a) per_bnd = 0 the periodicity is positiv
          (b) per_bnd = 1 the periodicity is negativ
        '''
        logger.info("negative periodicity")
        for i in range(n_copy):
            bi = np.append(bi,(-1)**(i)*bag[1][:-1])
            xi = np.append(xi,bag[0][-1]*(i)+bag[0][:-1])
    else:
        logger.info("positive periodicity")
        for i in range(n_copy):
            bi = np.append(bi,bag[1][:-1])
            xi = np.append(xi,bag[0][-1]*(i)+bag[0][:-1])
            

    Y = np.fft.fft(bi[0:-1])
    a = 2*np.abs(Y[pole_pairs])/(len(bi[0:-1]))
    time2 = time.time()
    logger.info("%s: %d pole pairs B amp %f Dauer der Berechnung %f",
                filename, pole_pairs, a, (time2-time1)*1000.0)
    

    time3 = time.time()
    N = 2**10  # The DFT is most efficient when N is a power of 2
    phi = np.linspace(0, 2*np.pi, N)
    nphi = int(round(360/((bag[0][-1] - bag[0][0])/(len(bag[0]) - 1)))) + 1
    ntiles = (nphi-1)//(len(bag[0]) - 1)
    phitab = np.linspace(0, 2*np.pi, nphi)
    if pmod % 2:
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

    br = np.interp(phi, phitab, bx)
    npoles = ntiles*pmod

    # compute DFT from induction
    Y = np.fft.fft(br)
    freq = np.fft.fftfreq(N, d=phi[1]-phi[0])
    
    # find the peak (amplitude of base harmonic)
    i = np.argmax(np.abs(Y[:N//2]))
    a = 2*np.abs(Y[i])/N
    T0 = np.abs(1/freq[i])
    time4 = time.time()
    logger.info("%s: %s poles B amp %f Dauer der Berechnung %f",
                filename, npoles, a, (time4-time3)*1000.0)
    
    alfa0 = np.angle(Y[i])
    alfa = bag[0]/180*np.pi
    
    logger.info("model angle: %f ",model_angle)
    logger.info("pole_pitch: %f ",pole_pitch)
    logger.info("number of copies: %d ",n_copy)
    logger.info("periodicity: %d ",per_bnd)
    logger.info("pole pairs: %d ",pole_pairs)    

    return dict(Bamp=a,
                phi0=alfa0,
                pos=bag[0].tolist(),
                B=bag[1].tolist(),
                nue=np.arange(0, 9*npoles).tolist(),
                B_nue=(2*np.abs(Y[:9*npoles])/N).tolist(),
                B_fft=(a*np.cos(2*np.pi*alfa/T0+alfa0)).tolist(),
                BX = bi.tolist(),
                XX = xi.tolist())
