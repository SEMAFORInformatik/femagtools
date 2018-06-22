# -*- coding: utf-8 -*-
"""
    femagtools.losscoeffs
    ~~~~~~~~~~~~~~~~~~~~~

    Fitting methods for loss coeffs



"""
import numpy as np
import scipy.optimize as so


def pfe_jordan(f, B, ch, fh, cw, fw, fb, fo, Bo):
    return (ch*(f/fo)**fh + cw*(f/fo)**fw)*(B/Bo)**fb


def pfe_steinmetz(f, B, cw, fw, fb, fo, Bo):
    return cw*(f/fo)**fw * (B/Bo)**fb


def fitsteinmetz(f, B, losses, Bo, fo):
    """fit coeffs of
    losses(f,B)=cw*(f/fo)**alfa*(B/Bo)**beta
    returns (cw, alfa, beta)
    """
    z = []
    for i, fx in enumerate(f):
        if fx:
            if isinstance(B[0], float):
                z += [(fx, bx, y)
                      for bx, y in zip(B, np.array(losses).T[i])
                      if isinstance(y, float)]
            else:
                z += [(fx, bx, y)
                      for bx, y in zip(B[i], np.array(losses).T[i])
                      if y]
                
    fbx = np.array(z).T[0:2]
    y = np.array(z).T[2]

    fitp, cov = so.curve_fit(
        lambda x, cw, alpha, beta: pfe_steinmetz(
            x[0], x[1], cw, alpha, beta, fo, Bo),
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
                z += [(fx, bx, y)
                      for bx, y in zip(B, np.array(losses).T[i])
                      if y]
            else:
                z += [(fx, bx, y)
                      for bx, y in zip(B[i], np.array(losses).T[i])
                      if y]

    fbx = np.array(z).T[0:2]
    y = np.array(z).T[2]
    fitp, cov = so.curve_fit(lambda x, cw, alpha, ch, beta, gamma: pfe_jordan(
        x[0], x[1], cw, alpha, ch, beta, gamma, fo, Bo),
                             fbx, y, (1.0, 2.0, 1.0, 1.0, 1.0))
    return fitp
