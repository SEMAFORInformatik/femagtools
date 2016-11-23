#!/usr/bin/env python
#
import unittest
import tempfile
import femagtools.poc
import os

sinPars = dict(
    pocType="Function",
    shape_current="sin",
    phi_voltage_winding=[0.0, 120.0, 240.0],
    skew_angle=30.0,
    num_winding=3,
    key_winding=[1, 2, 3],
    num_skew_steps=3)

funPars = dict(
    pocType="Function",
    shape_current="sin",
    func_steps=34,
    phi_voltage_winding=[-80.0, 40.0, 160.0],
    num_winding=3,
    key_winding=[1, 2, 3],
    func_current=[0.0,
                   11.25,
                   11.25,
                   33.75,
                   33.75,
                   56.25,
                   56.25,
                   78.75,
                   78.75,
                   101.25,
                   101.25,
                   123.75,
                   123.75,
                   146.25,
                   146.25,
                   168.25,
                   168.25,
                   191.25,
                   191.25,
                   213.75,
                   213.75,
                   236.25,
                   236.25,
                   258.75,
                   258.75,
                   281.25,
                   281.25,
                   303.75,
                   303.75,
                   326.25,
                   326.25,
                   348.75,
                   348.75,
                   360.0]
    ,func_phi= [0.0,
                0.0,
                0.4,
                0.4,
                1.1,
                1.1,
                1.66,
                1.66,
                2.0,
                2.0,
                1.66,
                1.66,
                1.1,
                1.1,
                0.4,
                0.4,
                0.0,
                0.0,
                -0.4,
                -0.4,
                -1.1,
                -1.1,
                -1.66,
                -1.66,
                -2.0,
                -2.0,
                -1.66,
                -1.66,
                -1.1,
                -1.1,
                -0.4,
                -0.4,
                0.0,
                0.0])

funcPars = dict(
        func_current=[0.0,0.0,1.0,1.0,0.0,0.0,-1.0,-1.0,0.0,0.0]
        ,func_phi=[0.0,30.0,30.0,150.0,150.0,210.0,210.0,330.0,330.0,360.0]
        ,harmonic_id=[1,3,5,7,9,11,13,15,17,19]
        ,func_steps=10 )


class PocFileTest(unittest.TestCase):

    def test_write_poc(self):
        poc = femagtools.poc.Poc(360, sinPars)
        filename = tempfile.mkstemp()[1]
        with open(filename, 'w') as f:
            poc.writefile(f)

        expected = '3\n1\n2\n3\n0.0\n120.0\n240.0\n360\nsin\n30.0\n3\n\n'
        with open(filename) as f:
            result = f.read()
        self.assertEqual(result, expected)

    def _createPoc(self, filename):
        testdir = os.path.split(__file__)[0]
        if not testdir:
            testdir = '.'
        return femagtools.poc.Poc(os.path.join(testdir, filename))

    def test_read_poc(self):
        poc = self._createPoc(filename='data/test.poc')
        expected = sinPars
        expected['pole_pitch'] = 360.0
        self.assertEqual(poc.getProps(), sinPars)

#    def test_read_poc_with_fun( self ):
#        poc=self._createPoc(filename='2p_sin.poc')
#        self.assertEqual(poc.getProps(), funPars)

if __name__ == '__main__':
    unittest.main()
