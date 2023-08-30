"""
    femagtools.vbf
    ~~~~~~~~~~~~~~

    Manage VBF magnetizing curve data files



"""
import sys
import femagtools.losscoeffs as lc
import logging

logger = logging.getLogger(__name__)

fo = 50.
Bo = 1.5


class Reader(object):

    def __init__(self, filename, filecontent=None):
        self.losses = {}
        # filecontent is used
        if filecontent:
            content = [l.strip() for l in filecontent.split('\n')]
        else:
            # only filename
            with open(filename) as f:
                content = [l.strip() for l in f.readlines()]

        if content:
            self.losses['name'] = content[0]
            self.losses['fo'], self.losses['Bo'] = [float(s)
                                                    for s in content[1].split()]
            # Ignore the next line
            self.losses['f'] = [float(s) for s in content[3].split()]
            self.losses['B'] = []
            self.losses['pfe'] = []
            for l in content[4:]:
                values = [float(s) for s in l.strip().split()]
                if len(values) > 1:
                    self.losses['B'].append(values[0])
                    self.losses['pfe'].append(
                        [v if v > 0 else None for v in values[1:]])
        logger.info("%s fmax %5.1f Bmax %3.2f", filename,
                    max(self.losses['f']), max(self.losses['B']))

    def __getitem__(self, index):
        return self.losses[index]
    
    def getLossValues(self):
        return self.losses


def read(filename, filecontent=None):
    """read VBF file and return dict of content"""
    vbf = Reader(filename, filecontent=filecontent)
    return vbf.getLossValues()
    

if __name__ == "__main__":
    import matplotlib.pylab as pl
    import femagtools.plot
    
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    losses = read(filename)
    print(losses)
    
    femagtools.plot.felosses(losses,
                             lc.fitjordan(
                                 losses['f'],
                                 losses['B'],
                                 losses['pfe'],
                                 losses['Bo'],
                                 losses['fo']),
                             title=filename, log=False)
    pl.show()
