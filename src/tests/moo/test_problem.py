#!/usr/bin/env python
#
import unittest
from femagtools import moo


class ProblemTest(unittest.TestCase):
    def test_compare_fc(self):
        prob = moo.Problem(0, 0, 2)
        f1 = (0, 0)
        f2 = (1, 1)
        c1 = []
        c2 = []
        self.assertTrue(prob.compare_fc(f1, c1, f2, c2))
        self.assertFalse(prob.compare_fc(f2, c2, f1, c1))
        self.assertFalse(prob.compare_fc(f2, c2, f2, c2))

if __name__ == '__main__':
    unittest.main()
