"""
    femagtools.tks
    ~~~~~~~~~~~~~~

    Manage TKS magnetizing curve data files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import sys
import numpy as np
import os
import re
import codecs
import femagtools.losscoeffs
import json
import logging

logger = logging.getLogger(__name__)

MUE0 = 4e-7*np.pi  # 1.2566371E-06
fo = 50.
Bo = 1.5
numPattern = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)')

MC1_MIMAX = 50
MC1_NIMAX = 50


def readlist(f):
    x = []
    while True:
        l = f.readline().strip()
        if not l:
            return list(zip(*x))
        vals = [float(n) for n in numPattern.findall(l.replace(',', '.'))]
        if vals:
            x.append(vals)


def pfe1(f, B, ch, fh, cw, fw, fb):
    return (ch*(f/fo)**fh + cw*(f/fo)**fw)*(B/Bo)**fb


def pfe2(f, B, cw, fw, fb):
    return cw*(f/fo)**fw * (B/Bo)**fb


class Reader:

    def __init__(self, filename):
        self.version_mc_curve = 0
        self.mc1_type = 1
        self.mc1_remz = 0.0
        self.mc1_recalc = 0
        self.mc1_bsat = 0.0
        self.mc1_bref = 0.0
        self.mc1_curves = 1

        self.mc1_fillfac = 1.0
        self.rho = 7.6
        self.MCURVES_MAX = 20
        self.MC1_NIMAX = 50
        self.MC1_MIMAX = 50
        self.mc1_db2 = [0. for I in range(self.MCURVES_MAX)]
#        self.mc1_energy = [[0 for I in range(self.MCURVES_MAX)]
#                           for K in range(self.MC1_MIMAX)]

        self.curve = []
        self.name = os.path.splitext(
            os.path.basename(filename))[0]
        self.mc1_title = 'ThyssenKrupp Steel'
        self.mc1_curves = 1
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
        self.fe_sat_mag = 2.15
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

        if self.losses and not np.isscalar(self.losses['B'][0]):
            import scipy.interpolate as ip
            z = femagtools.losscoeffs.fitjordan(
                self.losses['f'],
                self.losses['B'],
                self.losses['pfe'],
                self.Bo,
                self.fo)
            logger.info("Loss coeffs %s", z)
            self.ch = z[2]
            self.ch_freq = z[3]
            self.cw = z[0]
            self.cw_freq = z[1]
            self.b_coeff = z[4]
            
            z = femagtools.losscoeffs.fitsteinmetz(
                self.losses['f'],
                self.losses['B'],
                self.losses['pfe'],
                self.Bo,
                self.fo)
            
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

    def getValues(self):
            
        return {
            'name': self.name,
            'desc': self.mc1_title,
            'cversion': self.version_mc_curve,
            'ctype': self.mc1_type,
            'recalc': self.mc1_recalc,
            'remz': self.mc1_remz,
            'bsat': self.mc1_bsat,
            'bref': self.mc1_bref,
            'fillfac': self.mc1_fillfac,
            'curve': self.curve,
            'fo': self.fo,
            'Bo': self.Bo,
            'ch': self.ch,
            'ch_freq': self.ch_freq,
            'cw': self.cw,
            'cw_freq': self.cw_freq,
            'b_coeff': self.b_coeff,
            'rho': self.rho,
            'fe_sat_mag': self.fe_sat_mag,
            'losses': self.losses}
                                     
if __name__ == "__main__":
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()
            
    tks = Reader(filename)
    
    if tks.losses:
        import matplotlib.pylab as pl
        import numpy as np
        cw = tks.cw
        alpha = tks.cw_freq
        ch = tks.ch
        beta = tks.ch_freq
        gamma = tks.b_coeff
        
        for i, f in enumerate(tks.losses['f']):
            pfe = [p for p in np.array(tks.losses['pfe']).T[i] if p]
            pl.plot(tks.losses['B'], pfe1(f, np.array(tks.losses['B']),
                                          cw, alpha, ch, beta, gamma))
            pl.plot(tks. losses['B'][:len(pfe)], pfe,
                    marker='o', label="{} Hz".format(f))

        pl.title("Iron Losses " + filename)
        #pl.yscale('log')
        #pl.xscale('log')
        pl.xlabel("Induction [T]")
        pl.ylabel("Pfe [W/kg]")
        #pl.legend()
        pl.grid(True)
        #pl.savefig('tks.png')
        pl.show()

    mcv = tks.getValues()
    json.dump(mcv, sys.stdout)
