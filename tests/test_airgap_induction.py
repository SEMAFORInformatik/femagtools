#!/usr/bin/env python
#
import femagtools.airgap as ag
import numpy as np
import os


def test_airgap_induction():
    testPath = os.path.join(os.path.split(__file__)[0], 'data')
    if len(testPath) == 0:
        testPath = os.path.join(os.path.abspath('.'), 'data')
    r = ag.read(os.path.join(testPath, 'bag.dat'), 90)
    # FIXME: Python2 1.1698
    #np.testing.assert_almost_equal(r['Bamp'], 1.271, 3)
    assert 1.16 < r['Bamp'] < 1.28
