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
        if fx:
            if isinstance(B[0], float):
                z += [(fx, bx, y) for bx, y in zip(B, np.array(losses).T[i])
                      if isinstance(y, float)]
            else:
                z += [(fx, bx, y) for bx, y in zip(B[i], np.array(losses).T[i])
                      if y]
                
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
        if fx:
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
