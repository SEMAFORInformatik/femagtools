#!/usr/bin/env python
#
import femagtools.bch
import femagtools.machine
import femagtools.windings
import math
import pathlib
import pytest


@pytest.fixture
def data_dir():
    return pathlib.Path(__file__).with_name('data')


def test_char_ld_vector():
    m = dict(p=4,
             r1=0.0806, le=0.0, ls=0.0,
             wind_temp=20.0,
             ld=[0.0014522728, 0.0014522728],
             lq=[0.0038278836, 0.0032154],
             psim=[0.11171972, 0.11171972],
             i1=[80.0],
             beta=[-41.1, 0.0])

    pm = femagtools.machine.PmRelMachineLdq(3, m['p'],
                                            m['psim'],
                                            m['ld'],
                                            m['lq'],
                                            m['r1'],
                                            m['beta'],
                                            m['i1'])
    n = [50]
    u1 = 340.0
    T = [171.0]
    r = pm.characteristics(T, n, u1)

    assert r['i1'][0] == pytest.approx(79.873, rel=1e-1)
    assert r['u1'][0] == pytest.approx(304.090, rel=1e-1)
    assert r['beta'][0] == pytest.approx(-39.27, rel=1e-1)
    assert r['cosphi'][0] == pytest.approx(0.7584, rel=1e-1)


def test_char_ld_scalar():
    m = dict(p=4,
             r1=0.0806, le=0.0, ls=0.0,
             wind_temp=20.0,
             ld=0.0014522728,
             lq=0.0038278836,
             psim=0.11171972)

    pm = femagtools.machine.PmRelMachineLdq(3, m['p'],
                                            m['psim'],
                                            m['ld'],
                                            m['lq'],
                                            m['r1'])

    n = [50]
    u1 = 340
    T = [171]
    r = pm.characteristics(T, n, u1)

    assert r['i1'][0] == pytest.approx(79.2, rel=1e-1)
    assert r['u1'][0] == pytest.approx(321.77, rel=1e-1)
    assert r['beta'][0] == pytest.approx(-35.04, rel=1e-1)
    assert r['cosphi'][0] == pytest.approx(0.7225, rel=1e-1)


def test_char_psid_array():
    m = dict(p=6,
             r1=0.0, ls=0.0, le=0.0,
             wind_temp=50.0,
             u1=1500.0,
             rwdg=1.0,
             rlen=1.0,

             psid=[
                 [-2.2155, -1.8844, -1.5092, -1.10915, -0.69265,
                  -0.26285, 0.172865, 0.61355, 1.0584, 1.5064, 1.95335],
                 [-2.19485, -1.8669, -1.4938, -1.0941, -0.67725, -0.24738,
                     0.187985, 0.6272, 1.0682, 1.51025, 1.94215],
                 [-2.1329, -1.8102, -1.4462, -1.04965, -0.6349, -0.206885,
                     0.224, 0.6545, 1.0815, 1.49485, 1.8851],
                 [-2.03525, -1.7199, -1.36535, -0.98, -0.5754, -0.15967,
                     0.25879, 0.66955, 1.0626, 1.4385, 1.79655],
                 [-1.9173, -1.6093, -1.2698, -0.90125, -0.51485, -0.119,
                     0.275275, 0.65555, 1.01955, 1.3699, 1.701],
                 [-1.7885, -1.49205, -1.17005, -0.8246, -0.462, -0.09093,
                     0.27622, 0.62895, 0.9702, 1.2943, 1.6002],
                 [-1.66845, -1.38355, -1.077, -0.75285, -0.41545, -0.07042,
                     0.269535, 0.5971, 0.91455, 1.2166, 1.50255],
                 [-1.5582, -1.2824, -0.99505, -0.6902, -0.37485, -0.054705,
                     0.259315, 0.56385, 0.85995, 1.1445, 1.4154],
                 [-1.45355, -1.1977, -0.9233, -0.63665, -0.34132, -0.04326,
                     0.247065, 0.5327, 0.81027, 1.08395, 1.3391],
                 [-1.37165, -1.12385, -0.861, -0.59115, -0.313635, -0.032998,
                     0.23457, 0.50435, 0.76615, 1.0206, 1.27015],
                 [-1.2964, -1.0563, -0.8085, -0.55265, -0.290255, -0.0242655,
                     0.22365, 0.4781, 0.7273, 0.9695, 1.20785]],
             psiq=[
                 [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.67305, 0.71225, 0.7448, 0.7742, 0.8043, 0.83475, 0.86415,
                     0.89075, 0.9114, 0.92225, 0.91105],
                 [1.28765, 1.36675, 1.43815, 1.4997, 1.55505, 1.60545, 1.64955,
                     1.6842, 1.70135, 1.68245, 1.6212],
                 [1.8046, 1.91135, 2.00865, 2.093, 2.1637, 2.21795, 2.25435,
                     2.26765, 2.2484, 2.1966, 2.1182],
                 [2.2078, 2.32575, 2.43425, 2.5235, 2.59035, 2.6355, 2.6565,
                     2.64845, 2.61345, 2.55395, 2.46855],
                 [2.50495, 2.6215, 2.7293, 2.81575, 2.8784, 2.91515, 2.92425,
                     2.90745, 2.8665, 2.8007, 2.7104],
                 [2.7314, 2.8406, 2.9358, 3.01385, 3.0695, 3.10065, 3.1052,
                     3.08385, 3.03905, 2.9736, 2.8889],
                 [2.91375, 3.0086, 3.0926, 3.1605, 3.20775, 3.23295, 3.234,
                     3.2116, 3.1689, 3.10765, 3.03275],
                 [3.0632, 3.14755, 3.22035, 3.27915, 3.31905, 3.33935, 3.3376,
                     3.3159, 3.27845, 3.2298, 3.1556],
                 [3.19795, 3.2683, 3.3327, 3.38205, 3.4153, 3.4335, 3.42895,
                     3.409, 3.37505, 3.3271, 3.2662],
                 [3.31205, 3.38765, 3.43315, 3.47515, 3.5035, 3.521, 3.514,
                     3.4958, 3.46395, 3.4216, 3.36875]],
             Id=[-600.0, -540.0, -480.0, -420.0, -360.0,
                 -300.0, -240.0, -180.0, -120.0, -60.0, 0.0],
             iq=[0.0, 60.0, 120.0, 180.0, 240.0, 300.0,
                 360.0, 420.0, 480.0, 540.0, 600.0])

    pm = femagtools.machine.PmRelMachinePsidq(3, m['p'],
                                              m['psid'],
                                              m['psiq'],
                                              m['r1'],
                                              m['Id'], m['iq'])

    n = [16.67]
    u1 = 1400
    T = [5600.0]
    r = pm.characteristics(T, n, u1)

    assert r['i1'][0] == pytest.approx(207.74350, rel=1e-1)
    assert r['u1'][0] == pytest.approx(1192.88, rel=1e-1)
    assert r['beta'][0] == pytest.approx(-38.01, rel=1e-1)
    assert r['cosphi'][0] == pytest.approx(0.788, rel=1e-1)


def test_char_ldd_fieldweakening():
    m = dict(
        p=6,
        r1=0.0,
        ld=[0.005152],
        lq=[0.0104825],
        psim=[1.1298],
        beta=[0.0],
        i1=[0.0])

    pm = femagtools.machine.PmRelMachineLdq(3, m['p'],
                                            m['psim'],
                                            m['ld'],
                                            m['lq'],
                                            m['r1'],
                                            m['beta'],
                                            m['i1'])

    n = [16.67]
    u1 = 1400
    T = [5600.0]
    r = pm.characteristics(T, n, u1)

    assert r['i1'][0] == pytest.approx(211.971, rel=1e-1)
    assert r['u1'][0] == pytest.approx(1263.721, rel=1e-1)
    assert r['beta'][0] == pytest.approx(-30.0, rel=1e-1)
    assert r['cosphi'][0] == pytest.approx(0.729, rel=1e-1)


def test_i1beta_char():
    m = dict(
        p=4,
        r1=0.055,
        ld=[0.0012634272, 0.0012634272],
        lq=[0.003158568, 0.0027257272],
        psim=[0.11171972, 0.11171972],
        beta=[-35.0, 0.0],
        i1=[100.0])

    pm = femagtools.machine.PmRelMachineLdq(3, m['p'],
                                            m['psim'],
                                            m['ld'],
                                            m['lq'],
                                            m['r1'],
                                            m['beta'],
                                            m['i1'])

    n = [50]
    u1 = 340
    T = [171]
    r = pm.characteristics(T, n, u1)

    assert r['i1'][0] == pytest.approx(84.9758, rel=1e-1)
    assert r['u1'][0] == pytest.approx(278.48, rel=1e-1)
    assert r['beta'][0] == pytest.approx(-38.02, rel=1e-1)
    assert r['cosphi'][0] == pytest.approx(0.774, rel=1e-1)


def test_torque_iqd():
    m1 = femagtools.machine.PmRelMachineLdq(3, 4,
                                            psim=0.11172,
                                            ld=1.2667e-3,
                                            lq=3.1477e-3)
    beta = -35 * 3.1415/180
    i1 = 100
    iq, id = femagtools.machine.iqd(beta, i1)
    assert m1.torque_iqd(iq, id) == pytest.approx(215.87, rel=1e-1)

    beta = [-35, 0]
    i1 = [100]
    ld = [0.0012667696, 0.0012667696]
    lq = [0.0031477052, 0.0027165356]
    psim = [0.1117197, 0.1117197]

    m2 = femagtools.machine.PmRelMachineLdq(3, 4,
                                            psim=psim,
                                            ld=ld,
                                            lq=lq,
                                            beta=beta,
                                            i1=i1)
    assert m2.torque_iqd(iq, id) == pytest.approx(215.87, rel=1e-1)


def test_psidq_shortcircuit():
    psid = [[-0.51364332, -0.48331104, -0.44353648,
             -0.35671764, -0.15149428, 0.16361048]]
    psiq = [[0.0]*6]
    id = [-1000.0, -800.0, -600.0, -400.0, -200.0, 0.0]
    iq = [0]*6
    r1 = 0.1
    ls = 1e-3
    m1 = femagtools.machine.PmRelMachinePsidq(3, 4,
                                              psid=psid,
                                              psiq=psiq,
                                              id=id,
                                              iq=iq,
                                              r1=r1, ls=ls)
    w1 = 500
    iqs, ids = m1.iqd_uqd(w1, 0, 0)
    assert femagtools.machine.betai1(
        iqs, ids) == pytest.approx((-1.768, 42.48), rel=1e-1)


@pytest.mark.skip("temporarily ignored due to strange github action failures")
def test_ldq_create(data_dir):
    bch = femagtools.bch.read(str(data_dir / 'ldq.BATCH'))
    pm = femagtools.machine.create(bch, r1=0, ls=0)
    # with pytest.raises(ValueError) as exc_info:
    #    iqd = pm.iqd_torque(250)

    iqd = pm.iqd_torque(120)
    beta, i1 = femagtools.machine.betai1(*iqd)
    assert beta*180/math.pi == pytest.approx(-30, rel=1e-1)
    assert i1 == pytest.approx(62.9, rel=1e-1)


def test_ldq_losses(data_dir):
    bch = femagtools.bch.read(str(data_dir / 'ldq-losses.BATCH'))
    pm = femagtools.machine.create(bch, r1=0, ls=0)
    f1 = 200
    i1 = 120
    beta = -45/180*math.pi
    assert pm.betai1_plfe1(beta, i1, f1) == pytest.approx(835.3, rel=1e-1)
    assert pm.betai1_plfe2(beta, i1, f1) == pytest.approx(51.2, rel=1e-1)


def test_invpark():
    w1 = 314.15
    w1t = [w1*t/500.0 for t in range(6)]
    iq = 1
    id = 0
    ia, ib, ic = femagtools.machine.invpark(w1t, iq, id)
    assert ia[0] == pytest.approx(-1.0)
    assert ib[0] == pytest.approx(0.5)
    assert ic[0] == pytest.approx(0.5)
    assert ic[-1] == pytest.approx(-0.5, rel=1e-3)


def test_psidq_create():
    bch = femagtools.bch.Reader()
    bch.type = 'Fast Psid-Psiq-Identification'
    bch.machine = dict(p=6)
    bch.psidq = dict(
        psid=[[-2.4, -1.1376,  0.301322,  1.80345],
              [-2.08654, -0.978283,  0.385186,  1.65732],
              [-1.66857, -0.753761,  0.34352,  1.37547],
              [-1.42619, -0.581669,  0.277301,  1.15571]],
        psiq=[[-5.02114e-03, -1.91e-04, -7.67306e-05, -2.017e-04],
              [1.95562,  2.32212,  2.45695,  2.33965],
              [2.87904,  3.14081,  3.20475,  3.03648],
              [3.30830,  3.50945,  3.54963,  3.41783]],
        iq=[0., 200., 400., 600.],
        id=[-600., -400., -200., 0.])

    ls = 0
    r1 = 0
    pm = femagtools.machine.create(bch, r1, ls, lfe=1, wdg=0.6)
    T = 4003
    n = 10
    w1 = 2*math.pi*n*bch.machine['p']
    U = 338
    iqx, idx = pm.iqd_torque(T)
    assert idx == pytest.approx(-196.815433, rel=1e-1)
    assert iqx == pytest.approx(303.313471074, rel=1e-1)
    uq, ud = pm.uqd(w1, iqx, idx)
    assert ud == pytest.approx(-522.735358830, rel=1e-1)
    assert uq == pytest.approx(213.623501918, rel=1e-1)
    iqx, idx, tq = pm.iqd_torque_umax(T, w1, U)
    assert idx == pytest.approx(-297.986735702, rel=1e-1)
    assert iqx == pytest.approx(248.31913889352, rel=1e-1)


def test_psidq_losses(data_dir):
    bch = femagtools.bch.read(str(data_dir / 'psidq-losses.BATCH'))
    pm = femagtools.machine.create(bch, r1=0, ls=0)
    f1 = 200
    iq = 120
    id = -120
    assert pm.iqd_plfe1(iq, id, f1) == pytest.approx(887.4, rel=1e-1)
    assert pm.iqd_plfe2(iq, id, f1) == pytest.approx(50.4, rel=1e-1)


def test_wdg_resistance():
    from femagtools.machine.utils import wdg_resistance
    wdg = femagtools.windings.Winding(
        {'Q': 36, 'm': 3, 'p': 2, 'yd': 8, 'l': 2})
    n = 17  # wires / coil side
    g = 2  # parallel groups
    aw = 1.5083e-6  # wire cross section area m2
    da1 = 0.1172  # bore diameter m
    hs = 19.33e-3  # slot height m
    lfe = 0.1172  # length of stator lamination stack m
    sigma = 58e6  # conductivity 1/Ohm m (copper at 20°C)
    assert pytest.approx(0.3156, rel=1e-1) == wdg_resistance(
        wdg, n, g, aw, da1, hs, lfe, sigma)
