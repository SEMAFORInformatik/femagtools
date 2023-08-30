#!/usr/bin/env python
#
import pytest
import femagtools.machine.sizing


def test_im():
    P = 1.5e3
    fs = 12e3
    lda = 0.9
    p = 4
    f1 = 100
    n = f1/p
    udc = 550
    r = femagtools.machine.sizing.im(P, n, p, udc=udc,
                                     sigmas=fs, Ba=0.77,
                                     cos_phi=0.8, eta=0.8,
                                     lda=0.9)
    assert round(r['outer_diam'], 3) == 0.19
    assert r['stator']['num_slots'] == 36


def test_spm():
    P = 1.5e3
    fs = 12e3
    lda = 0.9
    p = 4
    f1 = 100
    n = f1/p
    udc = 550

    r = femagtools.machine.sizing.spm(P, n, p, udc=udc,
                                      Hc=700, sigmas=fs, brem=1.1, Ba=0.77,
                                      cos_phi=0.7, eta=0.8, demag=1.7,
                                      lda=0.9)
    assert round(r['outer_diam'], 3) == 0.19
    assert r['stator']['num_slots'] == 24


def test_eesm():
    P = 10e3
    speed = 4440/60
    p = 4
    udc = 600
    r = femagtools.machine.sizing.eesm(P, speed, p, udc=udc, Q1=36)
    assert r['rotor']['rot_hsm']
