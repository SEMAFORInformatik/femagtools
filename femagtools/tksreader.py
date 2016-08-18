"""
    femagtools.vbfreader
    ~~~~~~~~~~~~~~~~~~~~

    Manage TKS magnetizing curve data files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import sys
import math
import numpy as np
import os
import codecs
import femagtools.losscoeffs


MUE0 = 4e-7*math.pi  # 1.2566371E-06
fo = 50.
Bo = 1.5

MC1_MIMAX = 50
MC1_NIMAX = 50


def step(l):
    if len(l) < 2:
        return 0  
    return max([l[i+1]-l[i] for i in range(len(l)-1)])


def shift(a, b):
    """return offset of b within a"""
    for n in range(len(a)-1):
        if b[0]-(a[n+1]-a[n])/2-a[n] < 0:
            return n
    return len(a)


def findNotNone(l):
    """return lower and upper indexes of not none values in list"""
    for i in range(len(l)):
        if l[i]:
            break
    for j in range(len(l)-1, -1, -1):
        if l[j]:
            break
    return (i, j)


def pfe1(f, B, ch, fh, cw, fw, fb):
    return (ch*(f/fo)**fh + cw*(f/fo)**fw)*(B/Bo)**fb


def pfe2(f, B, cw, fw, fb):
    return cw*(f/fo)**fw * (B/Bo)**fb


class TksReader:

    def __init__(self, filename):
        bhc=False
        losses=False
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
        self.mc1_db2   = [ 0. for I in range(self.MCURVES_MAX)  ]
        self.mc1_energy=[ [ 0 for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_MIMAX) ]

        self.curve=[]
        self.name=os.path.splitext(os.path.basename(filename))[0]
        self.mc1_title='ThyssenKrupp Steel'
        self.mc1_curves=1
        self.curve=[{}]
        self.curve[0]['hi'] = []
        self.curve[0]['bi'] = []
        self.fo=fo
        self.Bo=Bo
        self.ch = 0.0
        self.ch_freq = 0.0
        self.cw = None
        self.cw_freq = None
        self.b_coeff =  None
        self.rho = 7.6
        self.fe_sat_mag = None
        self.losses=dict(f=[],B=[],pfe=[])

        k=0
        with codecs.open(filename, encoding='latin1') as f:
            for l in f:
                if l.startswith('H(A/m)	B(T)'):
                    bhc=True
                elif l.startswith('Comment'):
                    try:
                        self.mc1_title = l.split(':')[1].strip().encode('utf-8')
                    except UnicodeEncodeError:
                        self.mc1_title = l.split(':')[1].strip()
                        
                elif l.startswith('Mass Density'):
                    d=l.split(':')[-1].strip()
                    self.rho = float(d.split()[0].replace(',','.'))
                    if d.split()[-1]=='kg/m^3':
                        self.rho/=1e3
                elif l.startswith('f='):
                    self.losses['f'].append(float(l.split('=')[-1].split()[0].replace(',','.')))
                    self.losses['B'].append([])
                    self.losses['pfe'].append([])
                    
                elif l.startswith('B(T)'):
                    losses=True
                elif bhc:
                    r=l.split()
                    if len(r)>1:
                        h,b = [float(x.replace(',','.')) for x in r[:2]]
                        if h>0:
                            k=k+1
                            #print "{} {} {} {}".format( h, b, MUE0*h+j, (1-(MUE0*h+j)/b)*100)
                            self.curve[0]['hi'].append(h)
                            self.curve[0]['bi'].append(b)
                    else:
                        k=0
                        bhc=False

                elif losses:
                    r=l.split()
                    if len(r)==2:
                        self.losses['B'][-1].append(float(r[0].replace(',','.')))
                        self.losses['pfe'][-1].append(float(r[1].replace(',','.')))
                    else:
                        losses=False
                        fref=self.losses['f'][-1]
                                                        
            b=[]
            for x in self.losses['B']:
                b += x
            b.sort()
            th=step(b)/3
            b=[group.mean() for group in np.split(b, np.where(np.diff(b) > th)[0]+1)]
            for i in range(len(self.losses['pfe'])):
                n=shift(b,self.losses['B'][i])
                self.losses['pfe'][i]=n*[None]+self.losses['pfe'][i]
        
        if len(self.losses['pfe'])>0:
            maxlen=max([len(pfe) for pfe in self.losses['pfe']])
            pfe=[[] for i in range(maxlen)]
            for i in range(maxlen):
                    for p in self.losses['pfe']:
                        if i<len(p):
                                pfe[i].append(p[i])
                        else:
                            pfe[i].append(None)

            self.losses['pfe']=pfe
            self.losses['B']=b
        
    def getValues( self ):
        return {
            'name':self.name,
            'desc': self.mc1_title,
            'cversion': self.version_mc_curve,
            #                'ni': self.mc1_ni,
            #                'mi': self.mc1_mi,
            'ctype': self.mc1_type,
            'recalc': self.mc1_recalc,
            'db2': self.mc1_db2,
            'remz': self.mc1_remz,
            'bsat': self.mc1_bsat,
            'bref': self.mc1_bref,
            'fillfac': self.mc1_fillfac,
            ###                'curves': self.mc1_curves,
            'curve' : self.curve,
            'energy':self.mc1_energy,
            
            'fo':self.fo,
            'Bo': self.Bo,
            'ch': self.ch,
            'ch_freq': self.ch_freq,
            'cw':self.cw,
            'cw_freq':self.cw_freq,
            'b_coeff':self.b_coeff,
            'rho':self.rho,
            'fe_sat_mag':self.fe_sat_mag,
            'losses':self.losses }
                                     
if __name__ == "__main__":
    import matplotlib.pylab as pl
    if len(sys.argv)==2 :
        filename=sys.argv[1]
    else:
        filename=sys.stdin.readline().strip()
            
    tks = TksReader(filename)

    mcv=tks.getValues()

    for p in mcv['losses']['pfe']:
        print( p )
        
    cw, alfa, beta = femagtools.losscoeffs.fit(mcv['losses']['f'],
                                               mcv['losses']['B'],
                                               mcv['losses']['pfe'], 1.5, 50.)
    
    n=100
    B=pl.np.linspace(0.1,2,n)
    
    for i in range(len(mcv['losses']['f'])):
        f=mcv['losses']['f'][i]
        pfe=np.array(mcv['losses']['pfe']).T[i]
        j,k= findNotNone(pfe)
        pl.plot( B, pfe2(f,B,cw,alfa,beta) )
        pl.plot( mcv['losses']['B'][j:k], pfe[j:k], marker='o', label="f1={} Hz".format(f) ) 

    pl.title("Iron Losses " + filename)
    pl.yscale('log')
    pl.xscale('log')
    pl.xlabel("Induction [T]")
    pl.ylabel("Pfe [W/kg]")
    pl.grid( True )
    pl.show()

    #print tks.curve[0]['a']
    
    #json.dump( mcv, sys.stdout )
