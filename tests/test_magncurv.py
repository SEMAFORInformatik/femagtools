#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import unittest
import tempfile
import femagtools
import os


mcvPars =[dict(
    name="TKS_NO_20"
    ,desc="PowerCore NO 20 ;ThyssenKrupp Steel Eur"
    ,cversion=0
    ,ctype=1
    ,recalc=1
    ,curve=[{"bi":[0.0,0.09029103,0.17855673,0.26691416,0.35827994],
             "bi2":[11,22,33],
             "nuer":[111,222,333],
             "a":[6,5,4],
             "b":[9,8,7],
             "hi":[1,2,3,4,5]}]
    ,remz=0.0
    ,fillfac=0.92000002
    ,bsat=0.0
    ,bref=0.0
    ,ch=4.0
    ,ch_freq=5.0
    ,cw=3.0
    ,cw_freq=2.0
    ,fo=50.0
    ,Bo=1.5
    ,b_coeff=1.0
    ,rho=7.6500001
    ,fe_sat_mag=2.15)
]


class MagnetizingCurveTest(unittest.TestCase):
    maxDiff=None

    def test_findById( self ):
        mcv=femagtools.MagnetizingCurve(mcvPars)
        result=mcv.find( 'TKS_NO_20' )
        expected=mcvPars[0]['name']
        self.assertEqual(result, expected)

    def test_writeFile( self ):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath='.'
        mcvwriter = os.path.abspath(os.path.join( testPath, 'mcvwriter' ))
        mcv=femagtools.MagnetizingCurve(mcvPars)
        result=mcv.writefile( 'TKS_NO_20', '.', mcvwriter )
        self.assertEqual(result.split('.')[0], mcvPars[0]['name'])
        out=[]
        filename=result
        with open('dummy.mc') as f:
            out=f.readlines()
            
        self.assertEqual(len(out), 49)
        os.remove('dummy.mc')

if __name__ == '__main__':
  unittest.main()
