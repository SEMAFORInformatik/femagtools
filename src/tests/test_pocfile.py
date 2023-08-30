#!/usr/bin/env python
#
import femagtools.poc
import pathlib

sinPars = dict(
    shape_current="sin",
    phi_voltage_winding=[0.0, 120.0, 240.0],
    skew_angle=30.0,
    key_winding=[1, 2, 3],
    num_skew_steps=3)

funPars = dict(
    shape_current="sin",
    phi_voltage_winding=[-80.0, 40.0, 160.0],
    key_winding=['1', '2', '3'],
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
        ,harmonic_id=[1,3,5,7,9,11,13,15,17,19])


def test_default():
    poles = 4
    pole_pitch = 2*360/poles
    poc = femagtools.poc.Poc(pole_pitch)
    assert poc.filename() == f'sin_{poles}p.poc'
    expected = ['3','1','2','3','0.0','120.0','240.0',f'{pole_pitch}','sin','0.0','0','']
    assert poc.content() == expected


def test_write_poc(tmpdir):
    poc = femagtools.poc.Poc(360, sinPars)
    filename = tmpdir / 'sin.poc'

    poc.write(filename)
    expected = '3\n1\n2\n3\n0.0\n120.0\n240.0\n360\nsin\n30.0\n3\n'
    assert filename.read() == expected


def test_read_poc():
    datadir = pathlib.Path(__file__).resolve().parent.joinpath('data')
    poc = femagtools.poc.Poc(str(datadir / 'test.poc'))
    expected = sinPars
    expected['pole_pitch'] = 360.0
    expected['pocType'] = 'Function'
    assert poc.getProps() == expected

def test_read_func_poc():
    datadir = pathlib.Path(__file__).resolve().parent.joinpath('data')
    poc = femagtools.poc.Poc(str(datadir / 'test_func.poc'))
    assert poc.getProps()['pole_pitch'] == 18.0
    assert poc.getProps()['pocType'] == 'fun'

#    def test_read_poc_with_fun( self ):
#        poc=self._createPoc('2p_sin.poc')
#        self.assertEqual(poc.getProps(), funPars)

def test_rec(tmpdir):
    key_winding=['1','2', '3']
    phi_voltage_winding=[30, 150, 270]
    p = 2
    poc = femagtools.poc.Poc(360/p,
                             dict(key_winding=key_winding,
                                  phi_voltage_winding=phi_voltage_winding,
                                  shape_current='rec'))
    filename = tmpdir / 'rec.poc'
    poc.write(filename)
    expected = '3\n1\n2\n3\n30\n150\n270\n180.0\nrec\n0.0\n0\n'
    assert filename.read() == expected

def test_har(tmpdir):
    key_winding=['1','2', '3']
    phi_voltage_winding=[30, 150, 270]
    p = 2
    poc = femagtools.poc.Poc(360/p,
                             dict(key_winding=key_winding,
                                  phi_voltage_winding=phi_voltage_winding,
                                  func_current=[1, 0, 0.3, 0, 0.1],
                                  func_phi=[0, 0, 0, 0, 0],
                                  pocType='har'))
    filename = tmpdir / 'har.poc'
    poc.write(filename)
    expected = '3\n1\n2\n3\n30\n150\n270\n180.0\nhar\n5\n1, 0\n0, 0\n0.3, 0\n0, 0\n0.1, 0\n0.0\n0\n'
    assert filename.read() == expected

def test_hsp(tmpdir):
    poc = femagtools.poc.HspPoc(harm=[1,21],
                                amp=[1, 0.01],
                                phi=[0, 0])
    poc.pole_pitch = 360/2
    #poc = femagtools.poc.Poc(360/p,
    #                         dict(harmonic_id=[1, 21],
    #                              func_current=[1, 0.01],
    #                              func_phi=[0, 0],
    #                              pocType='hsp'))
    filename = tmpdir / 'hsp.poc'
    poc.write(filename)

    expected = '3\n1\n2\n3\n0.0\n120.0\n240.0\n180.0\nhsp\n2\n1, 1, 0\n21, 0.01, 0\n0.0\n0\n'
    assert filename.read() == expected

def test_fun(tmpdir):
    p = 2
    key_winding=['1','2', '3']
    phi_voltage_winding=[30, 150, 270]
    poc = femagtools.poc.Poc(360/p,
                             dict(key_winding=key_winding,
                                  phi_voltage_winding=phi_voltage_winding,
                                  func_current=[0, 1, -1, 0],
                                  func_phi=[0, 45, 135, 180],
                                  pocType='fun'))
    filename = tmpdir / 'fun.poc'
    poc.write(filename)
    expected = '3\n1\n2\n3\n30\n150\n270\n180.0\nfun\n4\n0, 0\n45, 1\n135, -1\n180, 0\n0, 0\n45, 1\n135, -1\n180, 0\n0, 0\n45, 1\n135, -1\n180, 0\n0.0\n0\n'
