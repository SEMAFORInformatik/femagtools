#!/usr/bin/env python
#
import os
import femagtools.tks


def test_read_tks():
    testPath = os.path.split(__file__)[0]
    if not testPath:
        testPath = '.'
    filename = "data/TKS-M400-65A.txt"
    tks = femagtools.tks.Reader('{0}/{1}'.format(testPath, filename))
    assert tks['name'] == 'TKS-M400-65A'
    assert len(tks['losses']['f']) == 4
    assert len(tks['losses']['B']) == 14
    assert len(tks['losses']['pfe']) == 4
    assert len(tks['losses']['pfe'][0]) == 14
    assert tks['losses']['f'][-1] == 400.0
    assert tks['losses']['B'][-1] == 1.9
    # [1.71, 3.55, 5.98, 9.21, 13.17, 17.95, 23.87, 30.84, 39.24, 48.9, 60.23, 73.1,, None, None]
    assert tks['losses']['pfe'][3][-3:] == [73.1, None, None]

