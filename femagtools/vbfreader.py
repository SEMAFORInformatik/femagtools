"""
    femagtools.vbfreader
    ~~~~~~~~~~~~~~~~~~~~

    Manage VBF magnetizing curve data files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import sys
import numpy as np


def findNotNone(l):
    """return lower and upper indexes of not none values in list"""
    for i in range(len(l)):
        if l[i]:
            break
    for j in range(len(l)-1, -1, -1):
        if l[j]:
            break
    return (i, j)

fo = 50.
Bo = 1.5


def pfe1(f, B, ch, fh, cw, fw, fb):
    return (ch*(f/fo)**fh + cw*(f/fo)**fw)*(B/Bo)**fb


def pfe2(f, B, cw, fw, fb):
    return cw*(f/fo)**fw * (B/Bo)**fb


def logpfe2(f, B, cw, fw, fb):
    return np.log10(cw) + fw*np.log10(f/fo) * fb*np.log10(B/Bo)


class VbfReader:

    def __init__(self, filename):
        self.vbf = {}
        with open(filename) as f:
            self.vbf['name'] = f.readline().strip()
            self.vbf['fo'],self.vbf['Bo'] = [float(s)
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

    def fitAllCoeffs(self):
        betacoeffs = []
        # fit beta: pfe = cw_b*(B/Bo)**beta
        for i, fref in enumerate(self.vbf['f']):
            if fref > 0:
                pfe = np.array(self.vbf['pfe']).T[i]
                j, k = findNotNone(pfe)
                if j <= k:
                    y = [np.log10(p) for p in pfe[j:k+1]]
                    x = [np.log10(b/self.vbf['Bo'])
                         for b in self.vbf['B'][j:k+1]]
                    A = np.vstack([x, np.ones(len(x))]).T
                    beta, w = np.linalg.lstsq(A, y)[0]
                    betacoeffs.append(beta)

        # fit alfa: pfe = cw_f*(f/fo)**alfa
        alfacoeffs=[]
        for i in range( len(self.vbf['B']) ):
            bref=self.vbf['B'][i]
            pfe=np.array(self.vbf['pfe'])[i]
            j,k= findNotNone(pfe)
            if j<=k:
                y=[ np.log10(p) for p in pfe[j:k+1] ]
                x=[ np.log10(f/self.vbf['fo']) for f in self.vbf['f'][j:k+1]]
                A=np.vstack([x, np.ones(len(x))]).T
                alfa,cw=np.linalg.lstsq(A, y)[0]
                if alfa>1.2 and alfa<1.8:
                    alfacoeffs.append(alfa)

        if len(self.vbf['f'])>1:
            alfa=np.average(alfacoeffs)
        else:
            alfa=1.3
        beta=np.average(betacoeffs)
        # fit cw: pfe = cw * (f/fo)**alfa * (B/Bo)**beta
        cw=[]
        for i in range( len(self.vbf['f']) ):
            f=self.vbf['f'][i]
            for k in range( len(self.vbf['B']) ):
                B=self.vbf['B'][k]
                if i<len(self.vbf['pfe'][k]):
                    pfe=self.vbf['pfe'][k][i]
                    if pfe:
                        a=(f/self.vbf['fo'])**alfa*(B/self.vbf['Bo'])**beta
                        cw.append(pfe/a)
        return( np.average(cw), alfa, beta)

if __name__ == "__main__":
    import matplotlib.pylab as pl
    if len(sys.argv)==2 :
        filename=sys.argv[1]
    else:
        filename=sys.stdin.readline().strip()
            
    vbf = VbfReader(filename)

    cw, alfa, beta = vbf.fitAllCoeffs()

    n=100
    B=pl.np.linspace(0.1,2,n)
    
    for i in range(len(vbf.vbf['f'])):
        f=vbf.vbf['f'][i]
        pfe=[p for p in np.array(vbf.vbf['pfe']).T[i] if p>0]
        pl.plot( B, pfe2(f,B,cw,alfa,beta) )
        pl.plot( vbf.vbf['B'][:len(pfe)], pfe, marker='o', label="f1={} Hz".format(f) ) 

    pl.title("Iron Losses " + filename)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlabel("Induction [T]")
    pl.ylabel("Pfe [W/kg]")
    pl.grid( True )
    pl.show()
