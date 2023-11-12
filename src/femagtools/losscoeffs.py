# -*- coding: utf-8 -*-
"""Fitting methods for loss coeffs


"""
import numpy as np
import scipy.optimize as so
import scipy.interpolate as ip


def pfe_bertotti0(f, B, ch, cw, ce):
    return (ch*f + cw*f**2)*B**2 + ce*f**1.5*B**1.5


def pfe_bertotti1(f, B, ch, alpha, cw, ce):
    return ch*f*B**alpha + cw*f**2*B**2 + ce*f**1.5*B**1.5


def pfe_jordan(f, B, ch, fh, cw, fw, fb, fo, Bo):
    return (ch*(f/fo)**fh + cw*(f/fo)**fw)*(B/Bo)**fb


def pfe_steinmetz(f, B, cw, fw, fb, fo, Bo):
    return cw*(f/fo)**fw * (B/Bo)**fb


def fitsteinmetz(f, B, losses, Bo, fo, alpha0=1.0):
    """fit coeffs of
    losses(f,B)=cw*(f/fo)**alfa*(B/Bo)**beta
    returns (cw, alfa, beta)
    """
    if np.isscalar(f):
        fbx = [[f]*len(losses), B]
        y = losses
        fitp, cov = so.curve_fit(
            lambda x, cw, beta: pfe_steinmetz(
                x[0], x[1], cw, alpha0, beta, fo, Bo),
            fbx, y, (1.0, 2.0))
        fitp = np.insert(fitp, 1, alpha0)

    else:
        pfe = losses
        z = []
        for i, fx in enumerate(f):
            if fx:
                if isinstance(B[0], float):
                    z += [(fx, bx, y)
                          for bx, y in zip(B, pfe[i])
                          if isinstance(y, float)]
                else:
                    z += [(fx, bx, y)
                          for bx, y in zip(B[i], pfe[i]) if y]

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
    pfe = losses
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


def fit_bertotti0(f, B, losses):
    """fit coeffs of
    losses(f,B)=(ch*f + ch*f**2)*B**2 + ce*f**1.5*B**1.5
    returns (ch, cw, ce)
    """
    pb = [ip.CubicSpline(bi, pi)
          for bi, pi in zip(B, losses)]
    i0 = 0
    if np.isclose(f[i0], 0):
        i0 = 1
    wy = [pb[i+i0](1.0)/f[i+i0]
          for i in range(len(f[i0:]))]
    csf = ip.CubicSpline(np.sqrt(f[i0:]), wy)
    xx = np.linspace(0, np.sqrt(f[-1]), 10)
    z = np.polyfit(xx, csf(xx), 2)
    ch0 = z[2]
    cw0 = z[0]
    ce0 = z[1]

    v = []
    for k in range(len(losses[i0])):
        for fx, bx, p in zip(f[i0:], B[i0:], losses[i0:]):
            if k < len(p):
                v.append((fx, bx[k], p[k]/fx))
            else:
                break

    def wbert(f, b, ch, cw, cx):
        return (ch + cw*f)*b**2 + cx*f**0.5*b**1.5

    z = np.array(v).T
    fbx = z[0:2]
    y = z[2]
    fitp, cov = so.curve_fit(
        lambda x, ch, cw, ce: wbert(
            x[0], x[1], ch, cw, ce),
        fbx, y, (ch0, cw0, ce0))
    return fitp

def fit_bertotti1(f, B, losses):
    """fit coeffs of
    losses(f,B)=ch*f*B**alpha + ch*f**2*B**2 + ce*f**1.5*B**1.5
    returns (ch, alpha, cw, ce)
    """
    v = []
    i0 = 0
    if np.isclose(f[0], 0):
        i0 = 1
    for k in range(len(losses[i0])):
        for fx, bx, p in zip(f[i0:], B[i0:], losses[i0:]):
            if k < len(p):
                v.append((fx, bx[k], p[k]/fx))
            else:
                break

    def wbert(f, b, ch, alpha, cw, cx):
        return ch*b**alpha + cw*f*b**2 + cx*f**0.5*b**1.5

    z = np.array(v).T
    fbx = z[0:2]
    y = z[2]

    fitp, cov = so.curve_fit(
            lambda x, ch, alpha, cw, cx: wbert(
                x[0], x[1], ch, alpha, cw, cx),
            fbx, y, (1e-3, 2, 1e-3, 1e-3))
    return fitp
