import femagtools
import importlib
import os
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')
                    
'''
  Use only pre-defined models as examples.
  Test if C:user\<windows-user>\femag exists.
  DXF and FE-calculations have separate examples.
'''

models = ['statorBG-magnetSector',
          'stator1-magnetIron3',
          'stator1-magnetIron4',
          'stator1-magnetIron5',
          'stator1-magnetIronV',
          'stator2-magnetSector',
          'stator1-spoke',
          'stator4-magnetSector',
          'statorRotor3-magnetIron',
          'statorRotor3-ipm-fml']

logger = logging.getLogger("fslcreator")
workdir = os.path.join(os.path.expanduser('~'), 'femag')
logger.info("Femagtools Version %s Working Dir %s",
            femagtools.__version__, workdir)
for m in models:
    mod = importlib.import_module(m)
    logger.info("--> %s <--", m)
    with open(os.path.join(workdir, m+'.fsl'), 'w') as f:
        f.write('\n'.join(getattr(mod, 'create_fsl')()))
