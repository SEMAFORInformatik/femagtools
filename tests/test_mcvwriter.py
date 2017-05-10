#!/usr/bin/env python3
#
import unittest
import os
import femagtools.mcv
import logging
import logging.config

def setup_logging():
    logfile = 'logging.json'
    if os.path.exists( logfile ):
        with open(logfile, 'rt') as f:
            config = json.load(f)
    else:
        config=dict(
            version = 1,
            disable_existing_loggers = False,
            formatters = {
                'simple': {'format':
                '%(name)-12s %(levelname)-8s %(message)s'}
            },
            handlers = {
                'console': {'class': 'logging.StreamHandler',
                'formatter': 'simple',
                "stream": "ext://sys.stderr"}
            },
            loggers = {
                'root': {'handlers': ['console'],
                         'level': logging.WARN}
            },
            root = {
                "level": "DEBUG",
                "handlers": ["console"]
            }
        )
    logging.config.dictConfig(config)



class McvReaderTest(unittest.TestCase):
    def test_read_mcv(self):
        testPath = os.path.split(__file__)[0]
        if not testPath:
            testPath = '.'
        filename = "data/TKS_NO_20.MCV"
        filename_out = "data/TKS_NO_20_out.MCV"
        reader = femagtools.mcv.Reader()
        reader.readMcv('{0}/{1}'.format(testPath, filename))
        r = reader.get_results()
        self.assertEqual(r['desc'],
                         u'PowerCore\xae NO 20 ;ThyssenKrupp Steel Eur')
        self.assertEqual(len(r['curve'][0]['bi']), 24)
        self.assertEqual(r['curve'][0]['bi'][0], 0.0)
        self.assertAlmostEqual(r['curve'][0]['bi'][-1], 1.836, 3)

        # test mcv writer
        import json
        writer = femagtools.mcv.Writer(r)
        # writer.setData(r)
        writeMcvFile = '{0}/{1}'.format(testPath, filename_out)
        writer.writeMcv(writeMcvFile)
        self.assertNotEqual(writer, None)

        # TEST
        reader2 = femagtools.mcv.Reader()
        reader2.readMcv(writeMcvFile)
        
        logger.debug("MC Ttile: %s", reader2.mc1_title)
        logger.debug("MC Version: %d", reader2.version_mc_curve)
        logger.debug("MC Type: %d", reader2.mc1_type)
        logger.debug("MC numCurves: %d", reader2.mc1_curves)
        logger.debug("MC FillFac: %f", reader2.mc1_fillfac)
        logger.debug("MC LIST [%s]", [reader.mc1_ni[0], reader.mc1_mi[0], reader.mc1_type, reader.mc1_recalc, reader.mc1_db2[0]])
        logger.debug("MC LIST [%s]", [reader2.mc1_ni[0], reader2.mc1_mi[0], reader2.mc1_type, reader2.mc1_recalc, reader2.mc1_db2[0]])
        
        for i in ['version_mc_curve', 'mc1_title', 'mc1_curves', 'mc1_title']:
            self.assertAlmostEqual(eval('reader.'+i), eval('reader2.'+i), 3)
        for i in ['mc1_ni[0]', 'mc1_mi[0]', 'mc1_type', 'mc1_recalc', 'mc1_db2[0]']:
            self.assertAlmostEqual(eval('reader.'+i), eval('reader2.'+i), 3)
        for i in ['mc1_remz', 'mc1_bsat', 'mc1_bref', 'mc1_fillfac']:
            self.assertAlmostEqual(eval('reader.'+i), eval('reader2.'+i), 3)
        for i in ['curve[0]["hi"]', 'curve[0]["bi"]', 'curve[0]["bi2"]', 'curve[0]["nuer"]']:
            self.assertAlmostEqual(eval('reader.'+i), eval('reader2.'+i), 3)
        for i in ['curve[0]["a"]', 'curve[0]["b"]']:
            self.assertAlmostEqual(eval('reader.'+i), eval('reader2.'+i), 3)
        for i in ['fo', 'Bo', 'ch', 'ch_freq', 'cw',
                  'cw_freq', 'b_coeff', 'rho', 'fe_sat_mag']:
            self.assertAlmostEqual(eval('reader.'+i), eval('reader2.'+i), 3)


setup_logging()
logger = logging.getLogger('test_read_mcv')
if __name__ == '__main__':
    unittest.main()
