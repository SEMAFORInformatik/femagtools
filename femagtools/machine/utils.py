"""
  femagtools.machine.utils

  auxiliary module
"""
import numpy as np
import numpy.linalg as la
import scipy.interpolate as ip
import logging

logger = logging.getLogger(__name__)


def K(d):
    """space phasor transformation matrix
    (Inverse Park Transformation) T-1 * dq
    arguments:
      d: rotation angle

    returns transformation matrix
    """
    return np.array((
        (-np.cos(d), np.sin(d)),
        (-np.cos(d-2*np.pi/3), np.sin(d-2*np.pi/3)),
        (-np.cos(d+2*np.pi/3), np.sin(d+2*np.pi/3))))


def T(d):
    """space phasor transformation matrix
    (Park Transformation) T * abc
    arguments:
      d: rotation angle

    returns transformation matrix
    """
    return np.array((
        (-np.cos(d), -np.cos(d-2*np.pi/3), -np.cos(d+2*np.pi/3)),
        (np.sin(d), np.sin(d-2*np.pi/3), np.sin(d+2*np.pi/3))))/3*2


def invpark(a, q, d):
    """ convert a dq vector to the abc reference frame
    (inverse park transformation)

    Args:
        a: rotation angle
        d: value in direct axis
        q: value in quadrature axis
    """
    if np.isscalar(a) and np.isscalar(q) and np.isscalar(d):
        return np.dot(K(a), (q, d))
    if np.isscalar(q) and np.isscalar(d):
        return np.array([K(x).dot((q, d)) for x in a]).T
    return np.array([K(x).dot((y, z)) for x, y, z in zip(a, d, q)]).T


KTH = 0.0039  # temperature coefficient of resistance
TREF = 20.0  # reference temperature of resistance
EPS = 1e-13


def xiskin(w, temp, zeta):
    return zeta*np.sqrt(abs(w)/(2*np.pi)/(50*(1+KTH*(temp-TREF))))


def kskinl(xi, nl):
    xi2 = 2*xi
    nl2 = nl*nl
    if np.isscalar(xi):
        if xi < EPS:
            k = 1
        else:
            k = (3 / (nl2*xi2)*(np.sinh(xi2) - np.sin(xi2)) /
                 (np.cosh(xi2)-np.cos(xi2)) +
                 ((nl2-1)/(nl2*xi)*(np.sinh(xi)+np.sin(xi)) /
                          (np.cosh(xi)+np.cos(xi))))
    else:
        xi2 = xi2[xi > EPS]
        k = np.ones(np.asarray(xi).shape)
        k[xi > 1e-12] = (3 / (nl2*xi2)*(np.sinh(xi2) - np.sin(xi2)) /
                         (np.cosh(xi2)-np.cos(xi2)) +
                         ((nl2-1)/(nl2*xi)*(np.sinh(xi)+np.sin(xi)) /
                          (np.cosh(xi)+np.cos(xi))))
    return k


def kskinr(xi, nl):
    xi2 = 2*xi
    nl2 = nl*nl
    return xi*((np.sinh(xi2)+np.sin(xi2))/(np.cosh(xi2)-np.cos(xi2))) + \
        ((nl2-1) / 3 * xi2*((np.sinh(xi)-np.sin(xi)) /
                            (np.cosh(xi)+np.cos(xi))))


def wdg_resistance(wdg, n, g, aw, da1, hs, lfe, sigma=56e6):
    """return winding resistance per phase in Ohm
    Arguments:
    wdg: (Winding) winding
    n: (int) number of wires per coil side
    g: (int) number of parallel coil groups
    lfe: length of stator lamination stack in m
    aw: wire cross section area m2
    sigma: (float) conductivity of wire material 1/Ohm m
    """
    # mean length of one turn
    lt = 2.8*(da1+hs)/2*wdg.yd*2*np.pi/wdg.Q + 16e-3 + 2*lfe
    return wdg.turns_per_phase(n, g)*lt/sigma/aw/g


def wdg_inductance(wdg, n, g, da1, lfe, ag):
    """return winding inductance per phase in H
    Arguments:
    wdg: (Winding) winding
    n: (int) number of wires per coil side
    g: (int) number of parallel coil groups
    da1: bore diameter in m
    lfe: length of stator lamination stack in m
    ag: length of airgap in m
    """
    return wdg.inductance(n, g, da1, lfe, ag)


def resistance(r0, w, temp, zeta, gam=0, nh=1):
    """return conductor resistance
    Arguments:
    r0: (float) dc resistance
    w: (float) current frequency in rad
    temp: (float) conductor temperature in deg Celsius
    zeta: (float) skin effect coefficient (penetration depth)
    gam: (float) constant coefficient (0..1)
    nh: (int) number of vertical conductors in slot"""
    xi = xiskin(w, temp, zeta)
    if np.isscalar(xi):
        if xi < 1e-12:
            k = 1
        else:
            k = (gam + kskinr(xi, nh)) / (1. + gam)
    else:
        k = np.ones(np.asarray(w).shape)
        k[xi > 1e-12] = (gam + kskinr(xi[xi > 1e-12], nh)) / (1. + gam)
    return r0*(1.+KTH*(temp - TREF))*k
# return r0*(1.+KTH*(temp - TREF))*(gam + kskinr(xi, nh)) / (1. + gam)


def leakinductance(l0, w, temp, zeta, nl=1, pl2v=0.5):
    """return conductor leakage inductance
    Arguments:
    r0: (float) dc resistance
    w: (float) current frequency in rad
    temp: (float) conductor temperature in deg Celsius
    zeta: (float) skin effect coefficient (penetration depth)
    nh: (int) number of vertical conductors in slot
    pl2v: (float) variable coefficient (0..1) """
    return l0*(1.0+pl2v*(kskinl(
        xiskin(w, temp, zeta), nl)-1))


def betai1(iq, id):
    """return beta and amplitude of dq currents"""
    return (np.arctan2(id, iq),
            la.norm((id, iq), axis=0)/np.sqrt(2.0))


def iqd(beta, i1):
    """return qd currents of beta and amplitude"""
    return np.sqrt(2.0)*i1*np.array([np.cos(beta),
                                     np.sin(beta)])


def puconv(dqpar, p, NR, UR, IR):
    """convert dqpar to per unit
    arguments:
    dqpar: dict from ld-iq or psid-psiq identification
    p: pole pairs
    NR: ref speed in 1/s
    UR: ref voltage per phase in V
    IR: ref current per phase in A
    """
    WR = 2*np.pi*NR*p
    PSIR = UR/WR
    SR = 3*UR*IR
    if 'beta' in dqpar:
        dqp = dict(beta=dqpar['beta'], losses=dict())
        dqp['i1'] = np.array(dqpar['i1'])/IR
    elif 'iq' in dqpar:
        dqp = dict(iq=np.array(dqpar['iq)'])/IR*np.sqrt(2), losses=dict())
        dqp['id'] = np.array(dqpar['id'])/IR*np.sqrt(2)
    else:
        raise ValueError('invalid dqpar')
    for k in 'psid', 'psiq':
        dqp[k] = np.array(dqpar[k])/PSIR
    if 'losses' in dqpar:
        for k in ('magnet', 'styoke_hyst', 'styoke_eddy',
                  'stteeth_hyst', 'stteeth_eddy', 'rotor_hyst', 'rotor_eddy'):
            dqp['losses'][k] = np.array(dqpar['losses'][k])/SR
        dqp['losses']['speed'] = p*dqpar['losses']['speed']/WR
        dqp['losses']['ef'] = dqpar['losses']['ef']
        dqp['losses']['hf'] = dqpar['losses']['hf']
    return dqp


def dqpar_interpol(xfit, dqpars, ipkey='temperature'):
    """return interpolated parameters at temperature or exc_current

    Arguments:
      xfit -- temperature or exc_current to fit dqpars
      dqpars -- list of dict with id, iq (or i1, beta), Psid and Psiq values
      ipkey -- key (string) to interpolate
    """
    # check current range
    ckeys = (('i1', 'beta'), ('id', 'iq'))
    dqtype = 0
    fpip = {k: dqpars[0][k] for k in ckeys[dqtype]}
    fpip['losses'] = dict()
    for k in ckeys[dqtype]:
        curr = np.array([f[k] for f in dqpars], dtype=object)
        shape = curr.shape
        if curr.shape != (len(dqpars), len(curr[0])):
            raise ValueError("current range conflict")
        curr = curr.astype(float)
        if not np.array([np.allclose(curr[0], c)
                         for c in curr[1:]]).all():
            raise ValueError("current range conflict")

    try:
        speed = np.array([d['losses']['speed'] for d in dqpars])
        if (np.max(speed) - np.min(speed))/np.mean(speed) > 1e-3:
            raise ValueError("losses: speed conflict")
    except KeyError:
        pass

    sorted_dqpars = sorted(dqpars, key=lambda d: d[ipkey])
    x = [f[ipkey] for f in sorted_dqpars]
    for k in ('psid', 'psiq'):
        m = np.array([f[k] for f in sorted_dqpars]).T
        if len(x) > 2:
            fpip[k] = np.array(
                [[ip.UnivariateSpline(x, y, k=2)(xfit)
                  for y in row] for row in m]).T
        else:
            fpip[k] = ip.interp1d(
                x, m, fill_value='extrapolate')(xfit).T
    try:
        for k in ('styoke_hyst', 'stteeth_hyst',
                  'styoke_eddy', 'stteeth_eddy',
                  'rotor_hyst', 'rotor_eddy',
                  'magnet'):
            m = np.array([f['losses'][k] for f in sorted_dqpars]).T
            if len(x) > 2:
                fpip['losses'][k] = np.array(
                    [[ip.UnivariateSpline(x, y, k=2)(xfit)
                      for y in row] for row in m]).T
            else:
                fpip['losses'][k] = ip.interp1d(
                    x, m, fill_value='extrapolate')(xfit).T
            fpip['losses']['speed'] = dqpars[0]['losses']['speed']
            for f in ('hf', 'ef'):
                if f in dqpars[0]['losses']:
                    fpip['losses'][f] = dqpars[0]['losses'][f]
    except KeyError:
        pass
    return x, fpip
