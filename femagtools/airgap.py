# -*- coding: utf-8 -*-
"""
    femagtools.airgap
    ~~~~~~~~~~~~~~~~~

    Read airgap dat file

    :copyright: 2017 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import numpy as np


def read(filename):
    """read dat file with columns (phi, Br, Bphi)
    returns samples, values, amplitude and phase of base harmonic

    Args:
      filename: the ame of the file to be processed
    """
    bag = np.loadtxt(filename).T
    if len(bag) < 3:
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
    if i == 0 or abs((br[-1] - br[0])/(np.max(br) - np.min(br))) > 1e-3:
        br = np.append(br, -1*br)
        Y = 2*np.fft.rfft(br)/len(br)
        freq = np.fft.fftfreq(len(br), d=dphi)
        i = np.argmax(abs(Y))

    if freq[i] > 0:
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

