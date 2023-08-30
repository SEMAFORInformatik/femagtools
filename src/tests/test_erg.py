#!/usr/bin/env python
#
import os
import femagtools.erg
import numpy as np


def read_erg(filename):
    testPath = os.path.join(os.path.split(__file__)[0], 'data')
    if len(testPath) == 0:
        testPath = os.path.join(os.path.abspath('.'), 'data')
    r = femagtools.erg.read(os.path.join(testPath, filename))
    return r


def test_read_erg():
    r = read_erg('ldlq.erg')
    assert len(r.keys()) == 14
    assert min(r['beta']) == -90.0
    assert max(r['beta']) == 0.0
    np.array(r['M_FE']).shape == (10, 10)

