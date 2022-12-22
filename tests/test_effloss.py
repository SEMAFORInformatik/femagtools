#!/usr/bin/env python
#
import pytest
import femagtools.machine.effloss


@pytest.fixture
def impars():
    return {'p': 2, 'm': 3, 'f1ref': 50, 'u1ref': 230.94, 'rotor_mass': 12.19, 'kfric_b': 1,
            'r1': 0.2, 'r2': 0.54459629551,
            'lsigma1': 0.0009996977709180593, 'lsigma2': 0.006836872295316844,
            'psiref': 0.728735006166921, 'wref': 314.1592653589793,
            'fec': 64.10000000000001, 'fee': 0, 'fexp': 7.0,
            'im': [2.043764851047857, 4.087529702095714, 6.1312945531435705,
                   8.175059404191428, 10.218824255239285],
            'psi': [0.27117452078113785, 0.5286659550273167, 0.7135568973434478,
                    0.8133415467307101, 0.8671946085010127]}


def test_imeffloss(impars):
    nmax = 8000/60
    T = 32.9
    u1 = 230
    temp = (120, 120)

    r = femagtools.machine.effloss.efficiency_losses_map(
        impars, u1, T, temp, nmax, npoints=(10, 4))
    assert r['T'] == pytest.approx(
        [-33.1, -0.6, 0.4, 32.7,
         -33.1, -0.6, 0.4, 32.7,
         -28.2, -0.6, 0.4, 25.9], abs=1e-1)
