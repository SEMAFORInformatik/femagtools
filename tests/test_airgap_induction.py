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


def test_airgap_induction2():
    testPath = os.path.join(os.path.split(__file__)[0], 'data')
    if len(testPath) == 0:
        testPath = os.path.join(os.path.abspath('.'), 'data')
    r = ag.read(os.path.join(testPath, 'bag2.dat'))
    np.testing.assert_almost_equal(r['Bamp'], 0.9272, 3)
    assert r['npoles'] == 32
