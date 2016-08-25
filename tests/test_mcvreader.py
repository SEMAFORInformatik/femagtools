#!/usr/bin/env python
#
import unittest
import os
import femagtools

class McvReaderTest(unittest.TestCase):
    def test_read_mcv( self ):
        testPath = os.path.split(__file__)[0]
        if not testPath: testPath='.'
        filename="data/TKS_NO_20.MCV"
        mcv = femagtools.McvReader()
        mcv.readMcv('{0}/{1}'.format(testPath, filename))
        r = mcv.get_results()
        self.assertEqual( r['desc'],
                          u'PowerCore\xae NO 20 ;ThyssenKrupp Steel Eur' ) 
        self.assertEqual( len(r['curve'][0]['bi']), 24 ) 
        self.assertEqual( r['curve'][0]['bi'][0], 0.0 ) 
        self.assertAlmostEqual( r['curve'][0]['bi'][-1], 1.836, 3 ) 

        
if __name__ == '__main__':
  unittest.main()
