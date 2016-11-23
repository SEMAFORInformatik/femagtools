#!/usr/bin/env python
#
import unittest
import os
import femagtools.vbf


class VbfReaderTest(unittest.TestCase):
    def test_read_vbf(self):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath = '.'
            
        filename = "data/m270_35.vbf"
        vbf = femagtools.vbf.Reader('{0}/{1}'.format(testPath, filename))
        r = vbf.getLossValues()
        self.assertEqual(r['name'], 'M235_35')
        self.assertEqual(len(r['f']), 7)
        self.assertEqual(len(r['B']), 19)
        self.assertEqual(len(r['pfe']), 19)
        self.assertEqual(len(r['pfe'][0]), 7)
        self.assertAlmostEqual(r['f'][-1], 2000.0, 1)
        self.assertAlmostEqual(r['B'][-1], 1.926, 3)
        self.assertAlmostEqual(r['pfe'][3][-1], 46.633, 3)

if __name__ == '__main__':
    unittest.main()
