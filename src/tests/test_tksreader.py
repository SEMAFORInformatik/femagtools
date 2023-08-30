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
    assert len(tks['losses']['B']) == 13
    assert len(tks['losses']['pfe']) == 4
    assert len(tks['losses']['pfe'][0]) == 13
    assert tks['losses']['f'] == [50.0, 100.0, 200.0, 400.0]
    assert round(tks['losses']['B'][-1], 1) == 1.8
    # 13.1279112, 17.92205587, 23.79869172, 30.84129721, 39.3434623, 48.97389844, 60.25629981, 72.88007909
    assert len([x for x in tks['losses']['pfe'][3] if x]) == 8
