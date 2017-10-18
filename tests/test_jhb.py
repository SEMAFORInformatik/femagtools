#!/usr/bin/env python
#
import unittest
import os
import femagtools.jhb


class JhbReaderTest(unittest.TestCase):
    def test_read_jhb(self):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath = '.'
        filename = "data/M270-50A_1000Hz_L.jhb"
        jhb = femagtools.jhb.Reader('{0}/{1}'.format(testPath, filename))
        self.assertEqual(jhb['name'], 'M270-50A_1000Hz_L')
        self.assertEqual(len(jhb['curve'][0]['hi']), 13)
        self.assertEqual(jhb['curve'][0]['bi'][-1], 1.17)

if __name__ == '__main__':
    unittest.main()
