#!/usr/bin/env python
#
import unittest
import os
import femagtools.mcv


class McvReaderTest(unittest.TestCase):
    def test_read_mcv(self):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath = '.'
        filename = "data/TKS_NO_20.MCV"
        reader = femagtools.mcv.Reader()
        reader.readMcv('{0}/{1}'.format(testPath, filename))
        r = reader.get_results()
        self.assertEqual(r['desc'],
                         u'PowerCore\xae NO 20 ;ThyssenKrupp Steel Eur')
        self.assertEqual(len(r['curve'][0]['bi']), 24)
        self.assertEqual(r['curve'][0]['bi'][0], 0.0)
        self.assertAlmostEqual(r['curve'][0]['bi'][-1], 1.836, 3)
        
    def test_read_oriented_mcv(self):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath = '.'
        filename = "data/V800-50A_aniso.MCV"
        reader = femagtools.mcv.Reader()
        reader.readMcv('{0}/{1}'.format(testPath, filename))
        r = reader.get_results()
        self.assertEqual(r['curve'][0]['angle'], 0)
        self.assertEqual(r['curve'][1]['angle'], 90)
        self.assertEqual(len(r['curve'][0]['bi']), 17)
        self.assertEqual(r['curve'][0]['bi'][0], 0.0)
        self.assertAlmostEqual(r['curve'][0]['bi'][-1], 1.973, 3)
        self.assertAlmostEqual(r['curve'][1]['bi'][-1], 1.9763, 3)
        self.assertAlmostEqual(r['curve'][0]['hi'][-1], 16073, 0)
        self.assertAlmostEqual(r['curve'][1]['hi'][-1], 18718, 0)

    def test_read_PM_mcv(self):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath = '.'
        filename = "data/FERRIT_20gC.MCV"
        mcv = femagtools.mcv.read(
            '{0}/{1}'.format(testPath, filename))
        self.assertEqual(mcv.mc1_type, 2)
        self.assertAlmostEqual(min(mcv.curve[0]['hi']), -560e3)
        self.assertAlmostEqual(max(mcv.curve[0]['hi']), -50e3)
        self.assertAlmostEqual(min(mcv.curve[0]['bi']), -1.10371697)
        self.assertAlmostEqual(max(mcv.curve[0]['bi']), 0.4)

        
if __name__ == '__main__':
    unittest.main()
