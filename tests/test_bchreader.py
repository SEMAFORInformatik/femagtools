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
        self.assertEqual(bch.version, '7.9.147')
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
        #self.assertAlmostEqual(bch.losses[0]['rotfe'], 0.000, 3)

        self.assertAlmostEqual(bch.lossPar['fo'][0], 50.0, 1)
        self.assertAlmostEqual(bch.lossPar['fo'][1], 50.0, 1)
        self.assertEqual(bch.get(('machine', 'p')), 2)
                                
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

        self.assertAlmostEqual(bch.lossPar['fo'][0], 50.0, 1)
        
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['stajo'],
                                       [[1, 100.0, 10.33, 15.804],
                                        [3, 300.0, 9.391, 142.234],
                                        [5, 500.0, 9.391, 395.094],
                                        [7, 700.0, 9.391, 774.383],
                                        [9, 900.0, 3.348, 455.591],
                                        [11, 1100.0, 2.971, 603.881],
                                        [13, 1300.0, 1.476, 419.063],
                                        [15, 1500.0, 0.882, 333.395]])
                               
        np.testing.assert_almost_equal(bch.losses[-1]['fft']['staza'],
                                       [[1, 100.0, 8.641, 13.065],
                                        [3, 300.0, 7.774, 117.587],
                                        [5, 500.0, 7.774, 326.631],
                                        [7, 700.0, 7.748, 637.999],
                                        [9, 900.0, 3.679, 500.663],
                                        [11, 1100.0, 2.915, 592.805],
                                        [13, 1300.0, 1.303, 370.023],
                                        [15, 1500.0, 0.626, 236.594]])
        
    def test_read_psidq(self):
        bch = self.read_bch('psidpsiq.BATCH')

        self.assertEqual(len(bch.torque_fft), 5)
        self.assertTrue('order' in bch.torque_fft[0])
        self.assertTrue('torque' in bch.torque_fft[0])
        self.assertEqual(len(bch.torque_fft[0]['torque']), 6)
        self.assertEqual(bch.torque_fft[0]['order'], [36, 48, 60, 72, 76, 84])

        self.assertEqual(sorted(bch.flux.keys()), ['1', '2', '3'])
        self.assertEqual(len(bch.flux['1']), 5)
        self.assertTrue('flux_k' in bch.flux['1'][0])
        self.assertEqual(len(bch.flux['1'][0]['flux_k']), 46)

        self.assertEqual(len(bch.torque), 5)
        self.assertEqual(len(bch.torque[-1]['torque']), 46)

        self.assertEqual(len(bch.psidq), 7)
        self.assertEqual(len(bch.psidq_ldq), 6)
        self.assertEqual(len(bch.psidq['psid']), 4)
        self.assertEqual(len(bch.psidq_ldq['ld']), 4)
        
        self.assertEqual(len(bch.psidq['losses']), 5)
        self.assertEqual(len(bch.psidq['losses']['styoke']), 4)

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

        self.assertEqual(len(bch.ldq['losses']), 5)
        self.assertEqual(len(bch.ldq['losses']['styoke']), 4)
        #self.assertTrue('i1' in bch.airgapInduction)
        #self.assertEqual(len(bch.airgapInduction['i1']), 5)
        #self.assertEqual(len(bch.airgapInduction['an']), 4)
        #self.assertEqual(len(bch.airgapInduction['an'][0]), 9)

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

    def test_dq(self):
        bch = self.read_bch('dq.BATCH')
        bch.get(['torque', 'torque']) == []
        bch.get(['linearForce[-1]', 'ripple_x']) == 0.0
        assert bch.get(['linearForce', 'ripple_z']) is None
        self.assertAlmostEqual(bch.dqPar['psid'][0], 2.7294321753800737, 5)
        self.assertAlmostEqual(bch.dqPar['psiq'][0], 1.0899999999999999, 5)
        
    def test_read_felosses(self):
        bch = self.read_bch('rel-felosses.BATCH')
        self.assertEqual(len(bch.losses), 4)
        self.assertEqual(bch.losses[-1]['stajo'], 4425.106)
        self.assertEqual(bch.losses[-1]['staza'], 7504.659)
    
if __name__ == '__main__':
    unittest.main()
