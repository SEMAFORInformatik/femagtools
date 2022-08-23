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
HBpattern = re.compile(r'H.+\s+B')
BPpattern = re.compile(r'B.+\s+P')


def readlist(section):
    x = []
    for l in section:
        if not l:
            break
        vals = [float(n) for n in numPattern.findall(l.replace(',', '.'))]
        if vals:
            x.append(vals)
    return list(zip(*x))


class Reader(object):

    def __init__(self, filename, filecontent=None):
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
        self.losses = dict(f=[], B=[])
        pfe = []

        # filecontent is used
        if filecontent:
            content = [l.strip() for l in filecontent.split('\n')]
        else:
            # only filename
            with codecs.open(filename, encoding='utf-8', errors='ignore') as f:
                content = [l.strip() for l in f.readlines()]
        if content:
            for i, l in enumerate(content):
                if HBpattern.match(l):
                    hbj = readlist(content[i+1:])
                    self.curve[0]['hi'] = hbj[0]
                    self.curve[0]['bi'] = hbj[1]

                elif l.startswith('Material Name'):
                    self.name = l.split(':')[1].strip()

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

                elif BPpattern.match(l):
                    b, p = readlist(content[i+1:])
                    self.losses['f'].append(fxref)
                    self.losses['B'].append(b)
                    pfe.append(p)
        logger.info("%s Bmax %3.2f", filename, max(self.curve[0]['bi']))

        if pfe and not np.isscalar(self.losses['B'][0]):
            colsize = max([len(p) for p in pfe])
            losses = [list(p) + [0]*(colsize-len(p)) for p in pfe]
            z = lc.fitjordan(
                self.losses['f'],
                self.losses['B'],
                losses,
                self.Bo,
                self.fo)
            logger.info("Jordan loss coeffs %s", z)
            self.ch = z[0]
            self.ch_freq = z[1]
            self.cw = z[2]
            self.cw_freq = z[3]
            self.b_coeff = z[4]

            z = lc.fitsteinmetz(
                self.losses['f'],
                self.losses['B'],
                losses,
                self.Bo,
                self.fo)
            logger.info("Steinmetz loss coeffs %s", z)

            self.losses['cw'] = z[0]
            self.losses['cw_freq'] = z[1]
            self.losses['b_coeff'] = z[2]

            self.losses['Bo'] = self.Bo
            self.losses['fo'] = self.fo

            # must normalize pfe matrix:
            (self.losses['B'],
             self.losses['pfe']) = femagtools.mcv.norm_pfe(self.losses['B'], pfe)

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


def read(filename, filecontent=None):
    """read Thyssen File TKS and return mc dict"""
    tks = Reader(filename, filecontent=filecontent)
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
        beta = tks.cw_freq
        ch = tks.ch
        alpha = tks.ch_freq
        gamma = tks.b_coeff

        femagtools.plot.felosses(tks.losses,
                                 (ch, alpha, cw, beta, gamma),
                                 #                                 (tks.losses['cw'],
                                 #                                  tks.losses['cw_freq'],
                                 #                                  tks.losses['b_coeff']),
                                 title=filename, log=False)
        pl.show()

    mcv = tks.getValues()
    json.dump(mcv, sys.stdout)
