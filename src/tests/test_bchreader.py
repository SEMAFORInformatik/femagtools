#!/usr/bin/env python
#
import unittest
import os
import femagtools.bch
from io import open
import numpy as np


class BchReaderTest(unittest.TestCase):

    def read_bch(self, filename):
        testPath = os.path.join(os.path.split(__file__)[0], 'data')
        if len(testPath) == 0:
            testPath = os.path.join(os.path.abspath('.'), 'data')
        r = femagtools.bch.Reader()
        with open('{0}/{1}'.format(testPath, filename),
                  encoding='latin1') as f:
            r.read(f)
        return r

    def test_read_cogging(self):
        bch = self.read_bch('cogging.BATCH')
        self.assertEqual(bch.version, '7.9.147 November 2012')
        self.assertEqual(bch.nodes, 2315)
        self.assertEqual(bch.elements, 3305)
        self.assertEqual(bch.quality, 100.0)

        self.assertEqual(len(bch.torque_fft), 1)
        self.assertEqual(len(bch.torque_fft[0]), 5)
        self.assertTrue('order' in bch.torque_fft[0])
        self.assertTrue('torque' in bch.torque_fft[0])
        self.assertEqual(len(bch.torque_fft[0]['torque']), 5)
        self.assertEqual(bch.torque_fft[0]['order'], [4, 12, 24, 36, 48])

        self.assertEqual(sorted(bch.flux.keys()), ['1', '2', '3'])
        self.assertEqual(sorted(bch.flux['1'][0].keys()),
                         sorted(['displ', 'voltage_four',
                                 'current_k', 'flux_k',
                                 'voltage_ir', 'displunit',
                                 'voltage_dpsi']))
        self.assertEqual(len(bch.flux['1'][0]['flux_k']), 61)
        self.assertEqual(bch.flux_fft['1'][0]['order'], [1, 3, 5, 7, 9, 11])

        self.assertEqual(len(bch.torque), 1)
        self.assertEqual(sorted(bch.torque[0].keys()),
                         sorted(['angle', 'force_y', 'force_x', 'torque',
                                 'current_1', 'ripple', 't_idpsi']))
        self.assertEqual(len(bch.torque[0]['torque']), 61)

        self.assertAlmostEqual(bch.losses[0]['winding'], 0.0, 1)
        self.assertAlmostEqual(bch.losses[0]['stajo'], 0.458, 2)
        self.assertAlmostEqual(bch.losses[0]['staza'], 0.344, 3)
        self.assertAlmostEqual(bch.losses[0]['magnetJ'], 0.006, 3)
        # self.assertAlmostEqual(bch.losses[0]['rotfe'], 0.000, 3)

        self.assertAlmostEqual(bch.lossPar['fo'][0], 50.0, 1)
        self.assertAlmostEqual(bch.lossPar['fo'][1], 50.0, 1)
        self.assertEqual(bch.get(('machine', 'p')), 2)
        np.testing.assert_almost_equal(bch.inertia, [0.230195e-3, 0.011774e-3])

    def test_read_sctest(self):
        bch = self.read_bch('sctest.BATCH')

        self.assertEqual(len(bch.torque_fft), 1)
        self.assertEqual(len(bch.scData['ia']), 134)
        self.assertAlmostEqual(bch.scData['ikd'], 0.0, 1)
        self.assertAlmostEqual(bch.scData['iks'], 1263.581, 2)
        self.assertAlmostEqual(bch.scData['tks'], 1469.736, 2)

    def test_read_pmsim(self):
        bch = self.read_bch('pmsim.BATCH')

        self.assertEqual(len(bch.torque_fft), 2)
        self.assertTrue('order' in bch.torque_fft[0])
        self.assertTrue('torque' in bch.torque_fft[0])
        self.assertEqual(len(bch.torque_fft[0]['torque']), 7)
        self.assertEqual(bch.torque_fft[1]['order'], [0, 12, 24, 30, 36, 42])

        self.assertEqual(sorted(bch.flux['1'][0].keys()),
                         sorted(['displ', 'voltage_four',
                                 'current_k', 'flux_k',
                                 'voltage_ir', 'displunit',
                                 'voltage_dpsi']))
        self.assertEqual(len(bch.flux['1'][0]['flux_k']), 46)

        self.assertEqual(len(bch.torque), 2)
        self.assertTrue('torque' in bch.torque[1])
        self.assertEqual(len(bch.torque[1]['torque']), 46)

        self.assertTrue('ld' in bch.dqPar)
        self.assertAlmostEqual(bch.dqPar['i1'][1], 49.992, 3)
        self.assertAlmostEqual(bch.dqPar['ld'][0], 9.9e-3, 6)
        self.assertAlmostEqual(bch.dqPar['ld'][0], 9.9e-3, 6)
        self.assertAlmostEqual(bch.dqPar['u1'][1], 358.38, 2)
        self.assertAlmostEqual(bch.dqPar['torque'][0], 65.3, 1)

        self.assertAlmostEqual(bch.machine['i1'], 50.0)

        self.assertAlmostEqual(bch.lossPar['fo'][0], 50.0, 1)

        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['order_el'],
                                       [1, 3, 5, 7, 9, 11, 13, 15])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['freq'],
                                       [100.0, 300.0, 500.0, 700.0, 900.0,
                                        1100.0, 1300.0, 1500.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['hyst'],
                                       [10.33, 9.391, 9.391, 9.391, 3.348,
                                        2.971, 1.476, 0.882])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['eddy'],
                                       [15.804, 142.234, 395.094, 774.383,
                                        455.591, 603.881, 419.063, 333.395])

        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['order_el'],
                                       [1, 3, 5, 7, 9, 11, 13, 15])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['freq'],
                                       [100.0, 300.0, 500.0, 700.0, 900.0, 1100.0, 1300.0, 1500.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['hyst'],
                                       [8.641, 7.774, 7.774, 7.748, 3.679, 2.915, 1.303, 0.626])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['eddy'],
                                       [13.065, 117.587, 326.631, 637.999, 500.663, 592.805, 370.023, 236.594])

    def test_read_pmsim_9(self):
        bch = self.read_bch('pmsim-9.BATCH')

        self.assertAlmostEqual(bch.machine['plfe'][0], 2540.2, 1)
        self.assertAlmostEqual(bch.machine['plfe'][1], 2020.5, 1)
        self.assertAlmostEqual(bch.dqPar['up'][0], 259.4, 1)

        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['order_mech'],
                                       [6, 18, 30, 42, 54, 90, 114])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['order_el'],
                                       [1.0, 3.0, 5.0, 7.0, 9.0, 15.0, 19.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['freq'],
                                       [400.0, 1200.0, 2000.0, 2800.0, 3600.0, 6000.0, 7600.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['hyst'],
                                       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['eddy'],
                                       [1637.884, 225.861, 93.969, 19.904, 6.661, 3.043, 1.752])

        assert [round(l*1e3, 4) for l in bch.dqPar['Lho']] == [0.5908, 0.6583]

    def test_read_relsim(self):
        bch = self.read_bch('relsim.BATCH')

        self.assertEqual(len(bch.torque), 1)
        self.assertTrue('torque' in bch.torque[0])
        self.assertAlmostEqual(np.mean(bch.torque[0]['torque']), 5.656, 2)
        self.assertAlmostEqual(bch.dqPar['u1'][1], 274.5, 1)
        self.assertAlmostEqual(bch.dqPar['torque'][0], 5.775, 1)

        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['freq'],
                                       [50.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['hyst'],
                                       [0.152])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['hyst'],
                                       [0.066])

    def test_read_pmsim_external(self):
        bch = self.read_bch('pmsim-external.BATCH')

        self.assertTrue('ld' in bch.dqPar)
        self.assertAlmostEqual(bch.dqPar['i1'][1], 49.992, 3)
        self.assertAlmostEqual(bch.dqPar['ld'][0], 0.86688e-3, 6)
        self.assertAlmostEqual(bch.dqPar['ld'][0], 0.86688e-3, 6)
        self.assertAlmostEqual(bch.dqPar['u1'][1], 2409.142, 2)
        self.assertAlmostEqual(bch.dqPar['torque'][0], 1137.92, 1)

        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['order_el'],
                                       [1, 3])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['freq'],
                                       [800.0, 2400.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['hyst'],
                                       [2619.555, 49.438])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo']['eddy'],
                                       [15512.529, 1186.523])

        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['order_el'],
                                       [1, 3, 5])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['freq'],
                                       [800.0, 2400.0, 4000.0])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['hyst'],
                                       [5688.175, 296.19, 0.989])
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza']['eddy'],
                                       [43864.352, 7108.561, 39.563])

    def test_read_psidq(self):
        bch = self.read_bch('psidpsiq.BATCH')

        self.assertEqual(len(bch.torque_fft), 10)
        self.assertTrue('order' in bch.torque_fft[0])
        self.assertTrue('torque' in bch.torque_fft[0])
        self.assertEqual(len(bch.torque_fft[0]['torque']), 7)
        self.assertEqual(bch.torque_fft[0]['order'],
                         [0, 4, 8, 12, 16, 20, 24])

        self.assertEqual(sorted(bch.flux.keys()), ['1', '2', '3'])
        self.assertEqual(len(bch.flux['1']), 10)
        self.assertTrue('flux_k' in bch.flux['1'][0])
        self.assertEqual(len(bch.flux['1'][0]['flux_k']), 16)

        self.assertEqual(len(bch.torque), 10)
        self.assertEqual(len(bch.torque[-1]['torque']), 16)

        self.assertEqual(len(bch.psidq), 7)
        self.assertEqual(len(bch.psidq_ldq), 6)
        self.assertEqual(len(bch.psidq['psid']), 3)
        self.assertEqual(len(bch.psidq_ldq['ld']), 3)

        self.assertEqual(len(bch.psidq['losses']), 17)
        self.assertEqual(len(bch.psidq['losses']['styoke']), 3)
        self.assertTrue('id' in bch.airgapInduction)
        self.assertEqual(bch.airgapInduction['id'],
                         [-200.0, -100.0, 0.0])
        self.assertEqual(len(bch.airgapInduction['Ba']), 3)
        self.assertEqual(len(bch.airgapInduction['Bm'][0]), 3)

    def test_read_ldq(self):
        bch = self.read_bch('ldq.BATCH')
        self.assertEqual(len(bch.torque_fft), 13)
        self.assertTrue('order' in bch.torque_fft[0])
        self.assertTrue('torque' in bch.torque_fft[0])
        self.assertEqual(len(bch.torque_fft[0]['torque']), 8)
        self.assertEqual(bch.torque_fft[0]['order'], [12, 36, 48, 56, 60,
                                                      72, 76, 84])

        self.assertEqual(sorted(bch.flux.keys()), ['1', '2', '3'])
        self.assertEqual(len(bch.flux['1']), 13)
        self.assertEqual(len(bch.flux['1'][0]), 7)
        self.assertTrue('flux_k' in bch.flux['1'][0])
        self.assertEqual(len(bch.flux['1'][0]['flux_k']), 46)

        self.assertEqual(len(bch.torque), 13)
        self.assertEqual(len(bch.torque[-1]['torque']), 46)

        self.assertEqual(len(bch.ldq['losses']), 7)
        self.assertEqual(len(bch.ldq['losses']['styoke']), 4)

        self.assertTrue('i1' in bch.airgapInduction)
        self.assertEqual(len(bch.airgapInduction['i1']), 3)
        self.assertEqual(len(bch.airgapInduction['an']), 4)
        self.assertEqual(len(bch.airgapInduction['an'][0]), 4)

    def test_read_pmsim2(self):
        bch = self.read_bch('PM_270_L8_001.BATCH')
        self.assertAlmostEqual(bch.dqPar['i1'][1], 70.0, 1)
        self.assertAlmostEqual(bch.dqPar['beta'][0], -38.0, 1)

    def test_read_linearforce(self):
        bch = self.read_bch('linearForce.BATCH')

        self.assertEqual(len(bch.linearForce), 1)
        self.assertEqual(len(bch.linearForce[0]['displ']), 26)
        self.assertEqual(bch.linearForce[0]['displ'][5], 15.0)
        self.assertEqual(bch.linearForce[0]['force_x'][7], -0.3439)
        self.assertEqual(bch.linearForce[0]['force_y'][2], 03107.0)
        self.assertEqual(bch.linearForce[0]['magnet_1'][13], 10.0)
        self.assertEqual(bch.linearForce_fft[0]['force'][0], 0.3483)
        self.assertEqual(bch.linearForce_fft[1]['force'][0], 3157.)

        self.assertEqual(len(bch.linearForce_fft), 2)
        self.assertEqual(len(bch.flux_fft), 3)

    def test_read_linmot_z(self):
        bch = self.read_bch('linmot_z.BATCH')
        self.assertEqual(len(bch.linearForce), 2)
        self.assertEqual(max(bch.linearForce[1]['force_z']), 4074.0)

    def test_dq(self):
        bch = self.read_bch('dq.BATCH')
        bch.get(['torque', 'torque']) == []
        bch.get(['linearForce[-1]', 'ripple_x']) == 0.0
        assert bch.get(['linearForce', 'ripple_z']) is None
#        self.assertAlmostEqual(bch.dqPar['psid'][0], 2.7294321753800737, 5)
#        self.assertAlmostEqual(bch.dqPar['psiq'][0], 1.0899999999999999, 5)
        self.assertAlmostEqual(bch.dqPar['psid'][0], 1.93, 5)
        self.assertAlmostEqual(bch.dqPar['psiq'][0], 0.77074639149333668, 5)

    def test_read_felosses(self):
        bch = self.read_bch('rel-felosses.BATCH')
        self.assertEqual(len(bch.losses), 4)
        self.assertEqual(bch.losses[-1]['stajo'], 4425.106)
        self.assertEqual(bch.losses[-1]['staza'], 7504.659)

    def test_read_pmsim_demag(self):
        bch = self.read_bch('PMREL-4p-skewed.BATCH')

        self.assertEqual(len(bch.demag), 9)
        self.assertEqual([-370.92, -2241.79, -2236.31],
                         [d['H_max'] for d in bch.demag if d['segment'] == 3])

    def test_read_characteristics(self):
        bch = self.read_bch('char.BATCH')
        self.assertEqual(len(bch.characteristics), 1)
        self.assertEqual(len(bch.characteristics[0].keys()), 19)
        self.assertEqual(len(bch.characteristics[0]['speed_torque']['n']), 16)

    def test_read_asterisks(self):
        bch = self.read_bch('PM-with-asterisks_001.BATCH')
        self.assertTrue(np.isnan(bch.nodes))
        self.assertAlmostEqual(bch.airgapInduction['an'][0][8][0], 0.0690, 1)
        self.assertAlmostEqual(bch.airgapInduction['an'][0][9][0], -0.9915, 1)

    def test_read_dist_leak(self):
        bch = self.read_bch('PM-4p-distleak.BATCH')
        self.assertTrue(bch.leak_dist_wind)
        self.assertEqual(bch.leak_dist_wind['nseg'], 4)


if __name__ == '__main__':
    unittest.main()
