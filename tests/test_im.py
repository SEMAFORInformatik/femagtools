import femagtools.machine.im
import pathlib
import pytest
import numpy as np


@pytest.fixture
def im():
    impars = {
        "p": 2, "m": 3, "f1type": 50, "u1type": 230, "rotor_mass": 6, "kfric_b": 1,
        "r1": 0.2, "r2": 0.5389, "lsigma1": 2e-3,
        "lsigma2": 4.4e-3, "kh": 10, "zeta2": 1,
        "fec": 120, "fee": 0, "fexp": 7.0,
        "iml": 5.317456519935536, "ims": 3.425169199057511, "mexp": 8.75033417091787}
    return femagtools.machine.im.InductionMachine(impars)


def test_im(im):
    u1 = 230
    f1 = 50
    s = 0.01
    torque = im.torqueu(2*np.pi*f1, u1, 2*np.pi*(1-s)*f1/im.p)
    assert pytest.approx(torque, rel=0.01) == 17.77
