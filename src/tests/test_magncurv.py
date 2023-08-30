#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import pytest
import femagtools.mcv
import os
import sys
import tempfile
import shutil


mcvPars = [dict(
    name="TKS_NO_20",
    desc="PowerCore NO 20 ;ThyssenKrupp Steel Eur",
    cversion=0,
    ctype=1,
    recalc=1,
    curve=[{"bi": [0.0, 0.09029103, 0.17855673, 0.26691416, 0.35827994],
            "bi2": [11, 22, 33],
            "nuer": [111, 222, 333],
            "a": [6, 5, 4],
            "b": [9, 8, 7],
            "hi": [1, 2, 3, 4, 5]}],
    remz=0.0,
    fillfac=0.92000002,
    bsat=0.0,
    bref=0.0,
    ch=4.0,
    ch_freq=5.0,
    cw=3.0,
    cw_freq=2.0,
    fo=50.0,
    Bo=1.5,
    b_coeff=1.0,
    rho=7.6500001,
    fe_sat_mag=2.15)
]


def test_findById():
    mcv = femagtools.mcv.MagnetizingCurve(mcvPars)
    result = mcv.find('TKS_NO_20')
    expected = mcvPars[0]['name']
    assert result == expected


def test_writeFile():
    dir = tempfile.mkdtemp()
    ext = '.MC' if sys.platform == 'win32' else '.MCV'
    mcv = femagtools.mcv.MagnetizingCurve(mcvPars)
    result = mcv.writefile('TKS_NO_20', dir)
    assert result == mcvPars[0]['name'] + ext

    mcv = femagtools.mcv.read(os.path.join(dir, result))
            
    assert pytest.approx(
        mcvPars[0]['curve'][0]['bi']) == mcv.get_results()['curve'][0]['bi']
    assert pytest.approx(
        mcvPars[0]['curve'][0]['hi']) == mcv.get_results()['curve'][0]['hi']

    shutil.rmtree(dir)


def test_writeFile_fillfac():
    ext = '.MC' if sys.platform == 'win32' else '.MCV'
    testPath = os.path.split(__file__)[0]
    mcv = femagtools.mcv.read(os.path.join(testPath, 'data/TKS_NO_20.MCV'))
    dir = tempfile.mkdtemp()
    m = femagtools.mcv.MagnetizingCurve(mcv)
    result = m.writefile('TKS_NO_20', dir, fillfac=0.5)
    assert result == mcvPars[0]['name'] + '-50' + ext
    
    bi = [0.0, 0.04908391088247299, 0.09705951809883118,
          0.14508341252803802, 0.1947421133518219,
          0.24495112895965576, 0.2946207821369171,
          0.3442765474319458, 0.39504921436309814,
          0.44526469707489014, 0.4969058334827423,
          0.5512984395027161, 0.5992029309272766,
          0.6502805352210999, 0.7150968909263611,
          0.754536509513855, 0.7817665934562683,
          0.8623831868171692, 0.9199547171592712,
          0.9600813984870911, 0.9857978820800781,
          1.0011695623397827, 1.0120362043380737,
          1.0208477973937988]

    filename = os.path.join(dir, result)
    mcv = femagtools.mcv.read(filename)
    
    assert pytest.approx(bi) == mcv.get_results()['curve'][0]['bi']
    shutil.rmtree(dir)


    
