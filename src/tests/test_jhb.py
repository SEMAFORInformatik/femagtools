#!/usr/bin/env python
#
import os
import femagtools.jhb


def test_read_jhb():
    testPath = os.path.split(__file__)[0]
    if not testPath:
        testPath = '.'
    filename = "data/M270-50A_1000Hz_L.jhb"
    jhb = femagtools.jhb.Reader('{0}/{1}'.format(testPath, filename))
    assert jhb['name'] == 'M270-50A_1000Hz_L'
    assert len(jhb['curve']) == 1
    assert len(jhb['curve'][0]['hi']) == 13
    assert jhb['curve'][0]['bi'][-1] == 1.17


def test_read_jhb_aniso():
    testPath = os.path.split(__file__)[0]
    if not testPath:
        testPath = '.'
    filename = "data/aniso.jhb"
    jhb = femagtools.jhb.Reader('{0}/{1}'.format(testPath, filename))
    assert jhb['name'] == 'aniso'
    assert len(jhb['curve']) == 2
    assert len(jhb['curve'][0]['hi']) == 39
    assert jhb['curve'][0]['bi'][-1] == 4.485

