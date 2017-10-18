#!/usr/bin/env python
#
import unittest
import os
import femagtools.tks


class TksReaderTest(unittest.TestCase):
    def test_read_tks(self):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath = '.'
        filename = "data/TKS-M400-65A.txt"
        tks = femagtools.tks.Reader('{0}/{1}'.format(testPath, filename))
        self.assertEqual(tks['name'], 'TKS-M400-65A')
        self.assertEqual(len(tks['losses']['f']), 4)
        self.assertEqual(len(tks['losses']['B']), 13)
        self.assertEqual(len(tks['losses']['pfe']), 13)
        self.assertEqual(len(tks['losses']['pfe'][0]), 4)
        self.assertAlmostEqual(tks['losses']['f'][-1], 400.0, 1)
        self.assertAlmostEqual(tks['losses']['B'][-1], 1.8, 3)
        self.assertAlmostEqual(tks['losses']['pfe'][3][-1], 30.841, 2)

if __name__ == '__main__':
    unittest.main()
