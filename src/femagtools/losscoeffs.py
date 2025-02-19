# -*- coding: utf-8 -*-
"""Fitting methods for loss coeffs


"""
import numpy as np
import scipy.optimize as so
import scipy.interpolate as ip
import logging


def pfe_bertotti(f, B, ch, alpha, cw, ce):
    return ch*f*B**alpha + cw*f**2*B**2 + ce*f**1.5*B**1.5


def pfe_jordan(f, B, ch, fh, cw, fw, fb, fo, Bo):
    return (ch*(f/fo)**fh + cw*(f/fo)**fw)*(B/Bo)**fb


def pfe_steinmetz(f, B, cw, fw, fb, fo, Bo):
    return cw*(f/fo)**fw * (B/Bo)**fb

def pfe_modified_steinmetz(f, B, ch, cw, alpha, beta, fo, Bo):
    return ch*(f/fo)*(B/Bo)**(alpha + B*beta) + cw*(f/fo)**2*(B/Bo)**2

def wbert(f, b, ch, cw, cx):
    return (ch + cw*f)*b**2 + cx*f**0.5*b**1.5

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

def fit_modified_steinmetz(f, B, losses, Bo, fo):
    """fit coeffs of modifeld steinmetz
    losses(f,B)=ch*(f/fo)*(B/Bo)**(alpha + B*beta) + cw*(f/fo)**2*(B/Bo)**2
    returns (ch, cw, alpha, beta)
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
    fitp, cov = so.curve_fit(lambda x, ch, cw, alpha, beta: pfe_modified_steinmetz(
        x[0], x[1], ch, cw, alpha, beta, fo, Bo),
        fbx, y, (1.0, 1.0, 1.0, 1.0))
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


def fit_bertotti(f, B, losses):
    """fit coeffs of
    losses(f,B)=(ch*f + ch*f**2)*B**2 + ce*f**1.5*B**1.5
    returns (ch, cw, ce)

    Args:
        f: list of n strictly increasing frequency values
        B: list of n or list of list of induction values (nxm)
        losses: list of list of fe loss values (nxm)
    """
    i0 = 0
    if np.isclose(f[i0], 0):
        i0 = 1
    v = []
    # collect losses and flux density along the frequency
    for k in range(len(losses[i0])):
        y = []
        bb = []
        if isinstance(B[0], tuple) or isinstance(B[0], list):
            for fx, bx, p in zip(f[i0:], B[i0:], losses[i0:]):
                if k < len(p):
                    y.append(p[k]/fx)
                    bb.append(bx[k])
        else:
            for fx, p in zip(f[i0:], losses[i0:]):
                if k < len(p):
                    y.append(p[k]/fx)
                    bb.append(B[k])
        j = len(y)
        if j > 2:
            v.append(np.array((f[i0:j+i0], bb, y)).T.tolist())

    z = np.array([b for a in v for b in a]).T
    fbx = z[0:2]
    y = z[2]
    fitp, _ = so.curve_fit(
        lambda x, ch, cw, ce: wbert(
            x[0], x[1], ch, cw, ce),
        fbx, y,
        bounds=[(0, 0, 0),
                (np.inf, np.inf, np.inf)])
    return fitp

def fit_bertotti0(f, B, losses, generate=False):
    """fit coeffs of
    losses(f,B)=(ch*f + ch*f**2)*B**2 + ce*f**1.5*B**1.5
    returns (ch, cw, ce)

    Args:
        f: list of n strictly increasing frequency values
        B: list of list of induction values (nxm)
        losses: list of list of fe loss values (nxm)
        generate: generates additional frequency samples if True
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
    df = 70
    # collect losses and flux density along the frequency
    for k in range(len(losses[i0])):
        y = []
        bb = []
        for fx, bx, p in zip(f[i0:], B[i0:], losses[i0:]):
            if k < len(p):
                y.append(p[k]/fx)
                bb.append(bx[k])
            else:
                break
        j = len(y)
        if j > 2:
            if generate:
                # generate additional samples to improve LM fit (experimental)
                nsteps = int(np.ceil((f[i0:][j-1] - f[i0])/df))
                fw = ip.CubicSpline(f[i0:j+1], y)
                bw = ip.CubicSpline(f[i0:j+1], bb)
                fx = np.linspace(f[i0], f[i0:][j-1], nsteps)
                v.append(np.array((fx, bw(fx), fw(fx))).T)
            else:
                v.append(np.array((f[i0:j+i0], bb, y)).T.tolist())


    def wbert(f, b, ch, cw, cx):
        return (ch + cw*f)*b**2 + cx*f**0.5*b**1.5


    z = np.array([b for a in v for b in a]).T
    fbx = z[0:2]
    y = z[2]
    guess = (ch0, cw0, ce0)
    fitp, _ = so.curve_fit(
        lambda x, ch, cw, ce: wbert(
            x[0], x[1], ch, cw, ce),
        fbx, y,
        bounds=[(0, 0, 0),
                (np.inf, np.inf, np.inf)])
    return fitp


def fit_bertotti1(f, B, losses):
    """fit coeffs of
    losses(f,B)=ch*f*B**alpha + ch*f**2*B**2 + ce*f**1.5*B**1.5
    returns (ch, alpha, cw, ce)
    """
    v = []
    i0 = 0
    if np.isclose(f[i0], 0):
        i0 = 1
    df = 70  # delta frequency
    # collect losses and flux density along the frequency
    for k in range(len(losses[i0])):
        y = []
        bb = []
        for fx, bx, p in zip(f[i0:], B[i0:], losses[i0:]):
            if k < len(p):
                y.append(p[k]/fx)
                bb.append(bx[k])
            else:
                break
        j = len(y)
        if j > 2:
            v.append(np.array((f[i0:j+1], bb, y)).T)

    def wbert(f, b, ch, alpha, cw, cx):
        return ch*b**alpha + cw*f*b**2 + cx*f**0.5*b**1.5

    z = np.array([b for a in v for b in a]).T
    fbx = z[0:2]
    y = z[2]

    fitp, cov = so.curve_fit(
            lambda x, ch, alpha, cw, cx: wbert(
                x[0], x[1], ch, alpha, cw, cx),
            fbx, y, (1e-3, 2, 1e-3, 1e-3))
    return fitp
