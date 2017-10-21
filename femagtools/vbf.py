"""
    femagtools.vbf
    ~~~~~~~~~~~~~~

    Manage VBF magnetizing curve data files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import sys
import numpy as np
import femagtools.losscoeffs as lc


fo = 50.
Bo = 1.5


class Reader(object):

    def __init__(self, filename):
        self.losses = {}
        with open(filename) as f:
            self.losses['name'] = f.readline().strip()
            self.losses['fo'], self.losses['Bo'] = [float(s)
                                                    for s in f.readline().strip().split()]
            # Ignore the next line
            f.readline()
            self.losses['f'] = [float(s) for s in f.readline().strip().split()]
            self.losses['B'] = []
            self.losses['pfe'] = []
            for l in f.readlines():
                values = [float(s) for s in l.strip().split()]
                if len(values) > 1:
                    self.losses['B'].append(values[0])
                    self.losses['pfe'].append(
                        [v if v > 0 else None for v in values[1:]])
            
    def __getitem__(self, index):
        return self.losses[index]
    
    def getLossValues(self):
        return self.losses


def read(filename):
    """read VBF file and return dict of content"""
    vbf = Reader(filename)
    return vbf.getLossValues()
    

if __name__ == "__main__":
    import matplotlib.pylab as pl
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    losses = read(filename)
    print(losses)
    
#    cw, alfa, beta = femagtools.losscoeffs.fitsteinmetz(
#        losses['f'],
#        losses['B'],
#        losses['pfe'],
#        losses['Bo'],
#        losses['fo'])

    B = pl.np.linspace(0.1, 2)
    
#    for i in range(len(losses['f'])):
#        f = losses['f'][i]
#        pfe = [p for p in np.array(losses['pfe']).T[i] if p]
#        pl.plot(B, pfe2(f, B, cw, alfa, beta))
#        pl.plot(vbf. vbf['B'][:len(pfe)], pfe,
#                marker='o', label="f1={} Hz".format(f))

#    pl.title("Iron Losses " + filename)
#    pl.yscale('log')
#    pl.xscale('log')
#    pl.xlabel("Induction [T]")
#    pl.ylabel("Pfe [W/kg]")
#    pl.grid(True)
#    pl.show()

    z = [lc.fitjordan(
        losses['f'],
        losses['B'],
        losses['pfe'],
        losses['Bo'],
        losses['fo']),
         lc.fitsteinmetz(
             losses['f'],
             losses['B'],
             losses['pfe'],
             losses['Bo'],
             losses['fo'])]
         

    cw, alpha, ch, beta, gamma = z[0]
#    cw, alpha, beta = z
#    print(cw, alpha, beta)
#    cw = 2.89
#    alpha = 1.43
#    beta = 1.85
    fo = losses['fo']
    Bo = losses['Bo']
    for i, f in enumerate(losses['f']):
        pfe = [p for p in np.array(losses['pfe']).T[i] if p]
        pl.plot(B, lc.pfe_jordan(f, B, cw, alpha, ch, beta, gamma, fo, Bo))
        pl.plot(losses['B'][:len(pfe)], pfe,
                marker='o', label="{} Hz".format(f))

    pl.title("Iron Losses " + filename)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlabel("Induction [T]")
    pl.ylabel("Pfe [W/kg]")
    pl.legend()
    pl.grid(True)
    pl.show()

    cw, alpha, beta = z[1]
    for i, f in enumerate(losses['f']):
        pfe = [p for p in np.array(losses['pfe']).T[i] if p]
        pl.plot(B, lc.pfe_steinmetz(f, B, cw, alpha, beta, fo, Bo))
        pl.plot(losses['B'][:len(pfe)], pfe,
                marker='o', label="{} Hz".format(f))

    pl.title("Iron Losses " + filename)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlabel("Induction [T]")
    pl.ylabel("Pfe [W/kg]")
    pl.legend()
    pl.grid(True)
    pl.show()
    
    for f in losses['f'][-3:]:
        pl.plot(B, lc.pfe_jordan(f, B, *z[0], fo=fo, Bo=Bo),
                label='{} Hz, Jordan'.format(f))
        pl.plot(B, lc.pfe_steinmetz(f, B, *z[1], fo=fo, Bo=Bo),
                label='{} Hz, Steinmetz'.format(f))
    pl.title("Iron Losses " + filename)
    pl.xlabel("Induction [T]")
    pl.ylabel("Pfe [W/kg]")
    pl.yscale('log')
    pl.xscale('log')
    pl.legend()
    pl.grid(True)
    pl.show()
