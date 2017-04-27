# -*- coding: utf-8 -*-
"""
    femagtools.losscoeffs
    ~~~~~~~~~~~~~~~~~~~~~

    Fitting methods for loss coeffs

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import numpy as np
import scipy.optimize as so


def fitsteinmetz(f, B, losses, Bo, fo):
    """fit coeffs of 
        losses(f,B)=cw*(f/fo)**alfa*(B/Bo)**beta
       returns (cw, alfa, beta)
    """
    z = []
    for i, fx in enumerate(f):
        z += [(fx, bx, y) for bx, y in zip(B, np.array(losses).T[i])
              if isinstance(y, float)]
    
    fbx = np.array(z).T[0:2]
    y = np.array(z).T[2]

    steinmetz = lambda x, cw, alpha, beta: (
        cw*(x[0]/fo)**alpha*(x[1]/Bo)**beta)

    fitp, cov = so.curve_fit(steinmetz,
                             fbx, y, (1.0, 1.0, 2.0))
    return fitp


def fitjordan(f, B, losses, Bo, fo):
    """fit coeffs of 
      losses(f,B)=(cw*(f/fo)**alpha + ch*(f/fo)**beta)*(B/Bo)**gamma
    returns (cw, alpha, ch, beta, gamma)
    """
    z = []
    for i, fx in enumerate(f):
        if isinstance(B[0], float):
            z += [(fx, bx, y) for bx, y in zip(B, np.array(losses).T[i])
                  if y]
        else:
            z += [(fx, bx, y) for bx, y in zip(B[i], np.array(losses).T[i])
                  if y]

    jordan = lambda x, cw, alpha, ch, beta, gamma: (
        (cw*(x[0]/fo)**alpha +
         ch*(x[0]/fo)**beta) *
        (x[1]/Bo)**gamma)

    fbx = np.array(z).T[0:2]
    y = np.array(z).T[2]
    fitp, cov = so.curve_fit(jordan,
                             fbx, y, (1.0, 2.0, 1.0, 1.0, 1.0))
    return fitp


if __name__ == "__main__":

    f = [100.0, 200.0, 400.0, 1000.0, 2000.0]
    pfe = [[0.15, 0.34, 0.67, 1.87, 4.37],
           [0.37, 0.76, 1.66, 4.72, 10.84],
           [0.62, 1.28, 2.78, 8.42, 19.27],
           [1.08, 2.24, 4.91, 14.9, 34.67],
           [2.03, 4.28, 9.32, 28.12, 69.33],
           [3.46, 7.29, 15.72, 47.91, 118.7]]
    B = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0]

    print(fitjordan(f, B, pfe, 1.5, 50.))
    print(fitsteinmetz(f, B, pfe, 1.5, 50.))
