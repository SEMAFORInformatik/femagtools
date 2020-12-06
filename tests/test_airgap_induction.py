#!/usr/bin/env python
#
import femagtools.airgap as ag
import numpy as np
import os


def test_airgap_induction():
    testPath = os.path.join(os.path.split(__file__)[0], 'data')
    if len(testPath) == 0:
        testPath = os.path.join(os.path.abspath('.'), 'data')
    r = ag.read(os.path.join(testPath, 'bag.dat'))
    np.testing.assert_almost_equal(r['Bamp'], 1.26914, 3)
    assert r['npoles'] == 8
    harms = [nue for nue, amp in zip(r['nue'],
                                     r['B_nue']) if amp > 1e-2]

    assert harms == [4, 12, 20, 28, 36, 44, 52, 60, 68]

def test_airgap_induction2():
    testPath = os.path.join(os.path.split(__file__)[0], 'data')
    if len(testPath) == 0:
        testPath = os.path.join(os.path.abspath('.'), 'data')
    r = ag.read(os.path.join(testPath, 'bag2.dat'))
    np.testing.assert_almost_equal(r['Bamp'], 0.9272, 3)
    assert r['npoles'] == 32
    harms = [nue for nue, amp in zip(r['nue'],
                                     r['B_nue']) if amp > 1e-2]

    assert harms == [16, 48, 80, 112, 144, 176, 208, 272]
