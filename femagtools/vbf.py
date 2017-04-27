"""
    femagtools.vbf
    ~~~~~~~~~~~~~~

    Manage VBF magnetizing curve data files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import sys
import numpy as np
import femagtools.losscoeffs


fo = 50.
Bo = 1.5


def pfe1(f, B, ch, fh, cw, fw, fb):
    return (ch*(f/fo)**fh + cw*(f/fo)**fw)*(B/Bo)**fb


def pfe2(f, B, cw, fw, fb):
    return cw*(f/fo)**fw * (B/Bo)**fb


def logpfe2(f, B, cw, fw, fb):
    return np.log10(cw) + fw*np.log10(f/fo) * fb*np.log10(B/Bo)


class Reader:

    def __init__(self, filename):
        self.vbf = {}
        with open(filename) as f:
            self.vbf['name'] = f.readline().strip()
            self.vbf['fo'], self.vbf['Bo'] = [float(s)
                                              for s in f.readline().strip().split()]
            # Ignore the next line
            f.readline()
            self.vbf['f'] = [float(s) for s in f.readline().strip().split()]
            self.vbf['B'] = []
            self.vbf['pfe'] = []
            for l in f.readlines():
                values = [float(s) for s in l.strip().split()]
                if len(values) > 1:
                    self.vbf['B'].append(values[0])
                    self.vbf['pfe'].append(
                        [v if v > 0 else None for v in values[1:]])
            
    def getLossValues(self):
        return self.vbf
    
if __name__ == "__main__":
    import matplotlib.pylab as pl
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()
            
    vbf = Reader(filename)

    print(vbf.getLossValues())
    
#    cw, alfa, beta = femagtools.losscoeffs.fitsteinmetz(
#        vbf.vbf['f'],
#        vbf.vbf['B'],
#        vbf.vbf['pfe'],
#        vbf.vbf['Bo'],
#        vbf.vbf['fo'])

    n = 100
    B = pl.np.linspace(0.1, 2, n)
    
#    for i in range(len(vbf.vbf['f'])):
#        f = vbf.vbf['f'][i]
#        pfe = [p for p in np.array(vbf.vbf['pfe']).T[i] if p]
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

    z = femagtools.losscoeffs.fitjordan2(
        vbf.vbf['f'],
        vbf.vbf['B'],
        vbf.vbf['pfe'],
        vbf.vbf['Bo'],
        vbf.vbf['fo'])

    cw, alpha, ch, beta, gamma = z
#    cw, alpha, beta = z
#    print(cw, alpha, beta)
#    cw = 2.89
#    alpha = 1.43
#    beta = 1.85
    for i, f in enumerate(vbf.vbf['f']):
        pfe = [p for p in np.array(vbf.vbf['pfe']).T[i] if p]
        pl.plot(B, pfe1(f, B, cw, alpha, ch, beta, gamma))
        pl.plot(vbf. vbf['B'][:len(pfe)], pfe,
                marker='o', label="{} Hz".format(f))

    pl.title("Iron Losses " + filename)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlabel("Induction [T]")
    pl.ylabel("Pfe [W/kg]")
    pl.legend()
    pl.grid(True)
    pl.show()
