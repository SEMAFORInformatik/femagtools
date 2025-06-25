import femagtools.machine.im
import pathlib
import pytest
import numpy as np

from femagtools.machine.im import InductionMachine

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

def test_characteristics():
    m = {'r1': 0.675074272807, 'ls': 0.0, 'le': 0.0, 'wdgcon': 1.0, 'bar_temp': 80.0, 
         'im': [0.815370706376, 3.15956148721, 5.50375226804, 7.84794304887, 10.1921338297], 
         'psi': [0.110531890489, 0.42241005873, 0.64470403922, 0.742820604784, 0.791593101916],  
         'lsigma1': 0.00282975049013, 'lsigma2': 0.00170976658221, 'r2': 0.526553117684, 
         'u1ref': 230.940107676, 'f1ref': 49.9998282765, 'rotor_mass': 10.0742734696, 'kfric_b': 1.0,
         'zeta2': 0.983497740554, 'pl2v': 0.801009432234, 'psiref': 0.717940875943, 'wref': 314.158186388, 
         'fec': 61.5, 'fee': 0.0, 'fexp': 7.0,  'kth1': 0.0039, 'kth2': 0.0039, 'p': 4, 'm': 3, 'u_1': 165.0}
    
    op = {'n': [0.0, 16.0, 41.0], 
          'T': [50.0, 47.5, 42.5]}

    tmech = True
    im = InductionMachine(m)
    u1 = m['u_1']
    m.pop('u_1')

    r = im.characteristics(op['T'], op['n'], u1,with_tmech=tmech)

    expected = {'losses': [237.01, 237.12, 236.46],
                  'eta': [0, 0.76, 0.88]}
    
    assert pytest.approx(expected['losses'], rel=1e-1) == r['losses']
    assert pytest.approx(expected['eta'], rel=1e-1) == r['eta']