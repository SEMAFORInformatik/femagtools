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
    pfe = np.asarray(losses).T
    z = []
    for i, fx in enumerate(f):
        if fx:
            if isinstance(B[0], float):
                z += [(fx, bx, y)
                      for bx, y in zip(B, pfe[i])
                      if isinstance(y, float)]
            else:
                z += [(fx, bx, y)
                      for bx, y in zip(B[i], pfe[i])
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
    losses(f,B)=(ch*(f/fo)**alpha + ch*(f/fo)**beta)*(B/Bo)**gamma
    returns (ch, alpha, cw, beta, gamma)
    """
    pfe = np.asarray(losses).T
    z = []
    for i, fx in enumerate(f):
        if fx:
            if isinstance(B[0], float):
                z += [(fx, bx, y)
                      for bx, y in zip(B, pfe[i])
                      if y]
            else:
                z += [(fx, bx, y)
                      for bx, y in zip(B[i], pfe[i])
                      if y]

    fbx = np.array(z).T[0:2]
    y = np.array(z).T[2]
    fitp, cov = so.curve_fit(lambda x, ch, alpha, cw, beta, gamma: pfe_jordan(
        x[0], x[1], ch, alpha, cw, beta, gamma, fo, Bo),
                             fbx, y, (1.0, 1.0, 1.0, 2.0, 1.0))
    return fitp
