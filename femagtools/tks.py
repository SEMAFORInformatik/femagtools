"""
    femagtools.tks
    ~~~~~~~~~~~~~~

    Manage TKS magnetizing curve data files



"""
import sys
import numpy as np
import os
import re
import codecs
import femagtools.losscoeffs as lc
import femagtools.mcv
import json
import logging

logger = logging.getLogger(__name__)

MUE0 = 4e-7*np.pi  # 1.2566371E-06
fo = 50.
Bo = 1.5
numPattern = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)')


def readlist(f):
    x = []
    while True:
        l = f.readline().strip()
        if not l:
            return list(zip(*x))
        vals = [float(n) for n in numPattern.findall(l.replace(',', '.'))]
        if vals:
            x.append(vals)


class Reader(object):

    def __init__(self, filename):
        self.version_mc_curve = 0
        self.mc1_type = femagtools.mcv.MAGCRV

        self.curve = []
        self.name = os.path.splitext(
            os.path.basename(filename))[0]
        self.mc1_title = 'ThyssenKrupp Steel'
        self.curve = [{}]
        self.curve[0]['hi'] = []
        self.curve[0]['bi'] = []
        self.fo = fo
        self.Bo = Bo
        self.ch = 0.0
        self.ch_freq = 0.0
        self.cw = None
        self.cw_freq = None
        self.b_coeff = None
        self.rho = 7.6
        self.losses = dict(f=[], B=[], pfe=[])

        with codecs.open(filename, encoding='utf-8', errors='ignore') as f:
            while True:
                l = f.readline().strip()
                if not l:
                    if not self.losses['B']:
                        self.losses = dict()
                    break
                
                if l.startswith('H(A/m)	B(T)') or l.startswith('H[A/m]	B[T]'):
                    h, b, j = readlist(f)
                    self.curve[0]['hi'] = h
                    self.curve[0]['bi'] = b
                    
                elif l.startswith('Comment'):
                    self.mc1_title = l.split(':')[1].strip()
                        
                elif l.startswith('Mass Density'):
                    d = numPattern.findall(l.replace(',', '.'))
                    self.rho = float(d[0])
                    if l.split()[-1] == 'kg/m^3':
                        self.rho /= 1e3
                        
                elif l.startswith('f='):
                    fref = numPattern.findall(l.replace(',', '.'))
                    fxref = float(fref[0])
                    b, p = readlist(f)
                    self.losses['f'].append(fxref)
                    self.losses['B'].append(b)
                    self.losses['pfe'].append(p)
                    
        logger.info("%s Bmax %3.2f", filename, max(self.curve[0]['bi']))

        if self.losses and not np.isscalar(self.losses['B'][0]):
            import scipy.interpolate as ip
            z = lc.fitjordan(
                self.losses['f'],
                self.losses['B'],
                self.losses['pfe'],
                self.Bo,
                self.fo)
            logger.info("Jordan loss coeffs %s", z)
            self.ch = z[2]
            self.ch_freq = z[3]
            self.cw = z[0]
            self.cw_freq = z[1]
            self.b_coeff = z[4]
            
            z = lc.fitsteinmetz(
                self.losses['f'],
                self.losses['B'],
                self.losses['pfe'],
                self.Bo,
                self.fo)
            logger.info("Steinmeth loss coeffs %s", z)
            
            self.losses['cw'] = z[0]
            self.losses['cw_freq'] = z[1]
            self.losses['b_coeff'] = z[2]
            self.losses['Bo'] = self.Bo
            self.losses['fo'] = self.fo
            
            # must normalize pfe matrix:
            bmin = max(list(zip(*(self.losses['B'])))[0])
            bmax = max([bx[-1] for bx in self.losses['B']])
            Bv = np.arange(np.ceil(10*bmin)/10.0,
                           (np.floor(10*bmax)+1)/10.0, 0.1)
            m = []
            for i, b in enumerate(self.losses['B']):
                pfunc = ip.interp1d(b, self.losses['pfe'][i], kind='cubic')
                bx = [x for x in Bv if x < b[-1]]
                m.append(pfunc(bx).tolist() + (len(Bv)-len(bx))*[None])
            self.losses['B'] = Bv.tolist()
            self.losses['pfe'] = list(zip(*m))

    def __getitem__(self, index):
        return self.__getattribute__(index)
    
    def getValues(self):
        """return values as mcv dict"""
        return {
            'name': self.name,
            'desc': self.mc1_title,
            'cversion': self.version_mc_curve,
            'ctype': self.mc1_type,
            'curve': self.curve,
            'fo': self.fo,
            'Bo': self.Bo,
            'ch': self.ch,
            'ch_freq': self.ch_freq,
            'cw': self.cw,
            'cw_freq': self.cw_freq,
            'b_coeff': self.b_coeff,
            'rho': self.rho,
            'losses': self.losses}


def read(filename):
    """read Thyssen File TKS and return mc dict"""
    tks = Reader(filename)
    return tks.getValues()

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()
            
    tks = Reader(filename)
    if tks.losses:
        import matplotlib.pylab as pl
        import femagtools.plot
        import numpy as np
        cw = tks.cw
        alpha = tks.cw_freq
        ch = tks.ch
        beta = tks.ch_freq
        gamma = tks.b_coeff

        femagtools.plot.felosses(tks.losses,
                                 (cw, alpha, ch, beta, gamma),
                                 title=filename, log=False)
        pl.show()
        
    mcv = tks.getValues()
    json.dump(mcv, sys.stdout)
