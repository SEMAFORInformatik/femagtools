# -*- coding: utf-8 -*-
"""
    femagtools.airgap
    ~~~~~~~~~~~~~~~~~

    Read airgap dat file



"""
import sys
import numpy as np
import logging

logger = logging.getLogger(__name__)


def read(filename):
    """read dat file with columns (phi, Br, Bphi)
    returns samples, values, amplitude and phase of base harmonic

    Args:
      filename: the ame of the file to be processed
    """
    bag = np.loadtxt(filename).T
    if len(bag) < 3:
        logger.warn("%s has incomplete content", filename)
        return(dict())

    phi = bag[0]
    br = bag[1]
    phisteps = [phi[i+1] - phi[i] for i in range(len(phi)-1)]
    dphi = sum(phisteps)/len(phisteps)
    freq = np.fft.fftfreq(len(br), d=dphi)

    # compute FFT from induction
    Y = 2*np.fft.rfft(br)/len(br)
    
    # find the peak and interpolate
    i = np.argmax(abs(Y))
    logger.debug("Len %d max phi %f max amp %f",
                 len(phi), phi[-1], abs(Y[i]))
    if i == 0 or abs((br[-1] - br[0])/(np.max(br) - np.min(br))) > 1e-2:
        br = np.append(br, -1*br)
        Y = 2*np.fft.rfft(br)/len(br)
        freq = np.fft.fftfreq(len(br), d=dphi)
        i = np.argmax(abs(Y))

    if freq[i] > 0:
        logger.info("%s: Period %f max phi %f B amp %f",
                    filename, 1/freq[i], len(br)*dphi, abs(Y[i]))
        a = abs(Y[i])
        T0 = abs(1/freq[i])
        x0 = np.angle(Y[i])
        phi = [ix*dphi for ix in range(len(br))]
        x = np.array([2*x/T0*np.pi for x in phi])

    return dict(Bamp=a,
                phi0=x0,
                pos=phi,
                B=br.tolist(),
                B_fft=(a*np.cos(x+x0)).tolist())

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()
            
    print(read(filename)['Bamp'])
