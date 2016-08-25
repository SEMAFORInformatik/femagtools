# -*- coding: utf-8 -*-
"""
    femagtools.magcurv
    ~~~~~~~~~~~~~~~~~~

    Creating and managing MCV/MC files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import json
import subprocess
import sys
import logging
import os.path
import numpy as np
import femagtools.losscoeffs

logger = logging.getLogger(__name__)

transl = dict(
    cversion='version_mc_curve',
    desc='mc1_title',
    
    recalc='mc1_recalc',
    remz='mc1_remz',
    ctype='mc1_type',
    rho='mc1_fe_spez_weigth',
    ch='mc1_ch_factor',
    cw='mc1_cw_factor',
    ch_freq='mc1_ch_freq_factor',
    cw_freq='mc1_cw_freq_factor',
    fillfac='mc1_fillfac',
    fillfac_old='mc1_fillfac_old',
    bref='mc1_bref',
    bsat='mc1_bsat',
    Bo='mc1_base_induction',
    b_coeff='mc1_induction_factor',
    fo='mc1_base_frequency',
    a='mc1_a',
    b='mc1_b',
    hi='mc1_hi',
    nuer='mc1_nuer',
    bi='mc1_bi',
    bi2='mc1_bi2',
    db2='mc1_db2',
    fe_sat_mag='mc1_fe_sat_magnetization'
    )
    
loss_transl = dict(
    f='frequency',
    B='induction',
    pfe='losses',
    cw='cw_m', 
    alfa='alfa_m',
    beta='beta_m', 
    fo='base_freq', 
    Bo='base_induct'
)

MC1_MIMAX = 50
MC1_NIMAX = 50
M_LOSS_INDUCT = 20
M_LOSS_FREQ = 20

MUE0 = 4e-7*np.pi  # 1.2566371E-06


def findNotNone(l):
    """return lower and upper indexes of not none values in list"""
    for i in range(len(l)):
        if l[i]:
            break
    for j in range(len(l)-1, -1, -1):
        if l[j]:
            break
    return (i, j)


class MagnetizingCurve(object):
    def __init__(self, mcvpar):
        """initialize this object either from
          a list of mcv objects or
          a single mcv or
          a directory"""
        self.mcv = {}
        if isinstance(mcvpar, list):
            logger.info("MagnetizingCurve is list")
            for m in mcvpar:
                if 'id' in m:
                    self.mcv[str(m['id'])] = m
                elif 'name' in m:
                    self.mcv[m['name']] = m
                    
        elif isinstance(mcvpar, dict):
            logger.info("MagnetizingCurve is dict")
            if 'id' in mcvpar:
                self.mcv[str(mcvpar['id'])] = mcvpar
            elif 'name' in mcvpar:
                self.mcv[mcvpar['name']] = mcvpar

        # Do not use unicode in PYTHON 3 as all strings are sequences of Unicode
        elif isinstance(mcvpar, str) or isinstance(mcvpar, unicode):
            self.mcdirectory = os.path.abspath(mcvpar)
            logger.info("MC Dir %s", self.mcdirectory)
        else:
            raise Exception("unsupported parameter type "+str(type(mcvpar)))
            
    def find(self, id):
        """find mcv by id or name"""
        try:
            return self.mcv[id]['name']
        except ValueError as ex:
            pass  # not found
        except KeyError as ex:
            try:
                ext = '.MC' if sys.platform == 'win32' else '.MCV'
                filename = ''.join((id, ext))
                logger.debug("search file %s in %s", filename,
                             self.mcdirectory)
                if os.access(os.path.join(self.mcdirectory,
                                          filename), os.R_OK):
                    return id
            except:
                pass
        logger.debug("search by name %s", id)
        m = self.find_by_name(id)
        return m['name'] if m else None
    
    def find_by_name(self, name):
        """find mcv by name"""
        for k in self.mcv.keys():
            if self.mcv[k]['name'] == name:
                    return self.mcv[k]
            # not found
        return None

    def recalc(self):
        for m in self.mcv:
            curve=self.mcv[m]['curve'][0]
            mi=MC1_MIMAX-2
            dh = curve['hi'][-1]-curve['hi'][-2]
            db = curve['bi'][-1]-curve['bi'][-2]
            dmue_d=db/dh
            dmue = curve['bi'][-1]/curve['hi'][-1]
            db = 3e-2*curve['bi'][-1]
            n3=1.5

            curve['muer']=[b/MUE0/h for b,h in zip( curve['bi'],curve['hi'] )]
            
            # extend bh-curve into saturation
            while dmue_d>1.01*MUE0 and dmue>1.5*MUE0:
                dmue_d=MUE0 + (dmue_d-MUE0)/np.sqrt(n3)
                curve['bi'].append(curve['bi'][-1]+db)
                curve['hi'].append(curve['hi'][-1]+db/dmue_d)
                curve['muer'].append(curve['bi'][-1]/MUE0/curve['hi'][-1])
                n3+=0.2
                dmue = curve['bi'][-1]/curve['hi'][-1]

            self.mcv[m]['db2'] = (curve['bi'][-1]**2 - curve['bi'][0]**2)/(mi-1)
            nuek0 = (curve['hi'][1]-curve['hi'][0])/\
              (curve['bi'][1]-curve['bi'][0])
            for j1 in range(len(curve['bi'])):
                    bk02=curve['bi'][j1]**2
                    if bk02>0: break
            curve['nuer']=[MUE0*nuek0]
            curve['bi2']=[bk02]
            curve['a']=[]
            curve['b']=[]

            bk1=0.0
            while bk1<=curve['bi'][-1]:
                    bk12 = bk02 + self.mcv[m]['db2']
                    bk1=np.sqrt(bk12)
                    j=2
                    while j<len(curve['bi']) and bk1 <= curve['bi'][j]:
                           j+=1
                    j-=1
                    bdel = curve['bi'][j] - curve['bi'][j1]
                    c1 = (curve['hi'][j] - curve['hi'][j1])/bdel
                    c2 = curve['hi'][j1] - c1*curve['bi'][j1]

                    nuek1 = c1 + c2/bk1

                    curve['a'].append(MUE0*(bk12*nuek0 - bk02*nuek1)/self.mcv[m]['db2'])
                    curve['b'].append(MUE0*(nuek1 - nuek0)/self.mcv[m]['db2'])
                    nuek0 = nuek1
                    bk02 = bk12

                    curve['nuer'].append(MUE0*nuek1)
                    curve['bi2'].append(bk12)

            curve['a'].append(1.0)
            curve['b'].append(MUE0*curve['hi'][-1]-curve['bi'][-1])
            
    def writefile( self, name, directory='.', writeproc='mcvwriter' ):
        """find magnetic curve by name or id and write binary file
        returns filename if found else None"""
        if not id:
            return None
        ext='.MC' if sys.platform=='win32' else '.MCV'
        mcv=self.find_by_name( name )
        if not mcv:
            filename=''.join((name,ext))
            try:
                import shutil
                logger.info( "Copy file %s", filename )
                shutil.copy( os.path.join(self.mcdirectory, filename), directory )
                return filename
            except:
                logger.error("MCV %s not found", str(filename) );
            return None

        filename=''.join((mcv['name'],ext))
        try:
            logger.info("create %s", str(filename) );
            proc=subprocess.Popen([writeproc], cwd=directory,
                                stdin=subprocess.PIPE,
                                stderr=subprocess.PIPE )

            proc.stdin.write( "&INPARAM\n".encode('latin1'))
            for k in mcv.keys():
                if k=='name':
                    proc.stdin.write( "filename='{}'\n".format(filename).encode('latin1'))
                elif k=='desc':
                    proc.stdin.write( "mc1_title='{}'\n".format(
                        mcv['desc']).encode('latin1') )
                elif k=='curve':
                    proc.stdin.write('mc1_curves={}\n'.format(len(mcv[k])).encode('latin1'))
                    for c in mcv['curve']:
                        for n in c.keys():
                            if not n in transl: continue
                            proc.stdin.write( "{}=".format(transl[n]).encode('latin1'))
                            if type(c[n]) is list: proc.stdin.write(
                                    ','.join(map(str,c[n])).encode('latin1') )
                            else: proc.stdin.write( str(c[n]).encode('latin1') )
                            proc.stdin.write( "\n".encode('latin1') )

                    if 'energy' in mcv['curve'][0] and len(mcv['curve'][0]['energy'])>0:
                        proc.stdin.write( "mc1_energy=" )
                        proc.stdin.write( ','.join(map(str,mcv['curve'][0]['energy'])))
                        proc.stdin.write( "\n" )

                elif k=='losses':
                    # find start index
                    bstart=0
                    for a in mcv[k]['pfe']:
                        for i in range(len(a)):
                            if a[i]:
                                if i>bstart:
                                    bstart=i
                                break
                    proc.stdin.write('N_freq={}\n'.format(len(mcv[k]['f'])).encode('latin1'))
                    proc.stdin.write('N_J_ind={}\n'.format(len(mcv[k]['B'])-bstart).encode('latin1'))
                    proc.stdin.write( "frequency=".encode('latin1'))
                    proc.stdin.write( ",".join(map(str,mcv[k]['f'])).encode('latin1') )
                    proc.stdin.write( "\n" )
                    proc.stdin.write( "induction=" )
                    proc.stdin.write( ",".join(map(str, mcv[k]['B'][bstart:]) ).encode('latin1') )
                    proc.stdin.write( "\n".encode('latin1') )
                    proc.stdin.write( "losses=".encode('latin1') )
                    pfeT=[]
#                    flen=10
#                    blen=20
                    cw=mcv[k]['cw']
                    alfa=mcv[k]['alfa']
                    beta=mcv[k]['beta']
                    fo=mcv[k]['fo']
                    Bo=mcv[k]['Bo']

                    pfe=[]
                    lower=0
                    # must replace all None by approx losses
                    for i in range( len(mcv['losses']['f']) ):
                        f=mcv['losses']['f'][i]
                        if f>0:
                            pfei=[p[i] if i<len(p) else None for p in mcv['losses']['pfe']]
                            m,n= findNotNone(pfei)
                            if m>lower: lower=m
                            if m<=n:
                                y=[ np.log10(p) for p in pfei[m:n+1] ]
                                x=[ np.log10(b/Bo) for b in mcv['losses']['B'][m:n+1]]
                                A=np.vstack([x, np.ones(len(x))]).T
                                beta,cw=np.linalg.lstsq(A, y)[0]
                                for j in range(n+1,len(pfei)):
                                    pfei[j]=10**cw*(mcv['losses']['B'][j]/Bo)**beta

                                pfe.append(pfei)
                    pfe+=[[0]*M_LOSS_INDUCT]*(M_LOSS_FREQ-len(mcv['losses']['f']))

                    for r in pfe:
                        a=list(r[lower:]) + [0]*(M_LOSS_INDUCT-len(r[lower:]))
                        #must copy last and cut "None" (null)
                        for i in range(len(a)-1,0,-1):
                            if a[i]>0: break
                            if a[i-1]>0:
                                a[i]=a[i-1]
                                break
                        pfeT+=a
                        
                    proc.stdin.write(','.join(map(str,pfeT)).encode('latin1'))
                    proc.stdin.write( "\n".encode('latin1') )
                    proc.stdin.write( "cw_m={}\n".format(mcv[k]['cw']).encode('latin1'))
                    proc.stdin.write( "alfa_m={}\n".format(mcv[k]['alfa']).encode('latin1'))
                    proc.stdin.write( "beta_m={}\n".format(mcv[k]['beta']).encode('latin1'))
                    proc.stdin.write( "base_freq={}\n".format(mcv[k]['fo']).encode('latin1'))
                    proc.stdin.write( "base_induct={}\n".format(mcv[k]['Bo']).encode('latin1'))
                    logger.info("has losses N_freq %d, N_J_ind %d", len(mcv[k]['f']), len(mcv[k]['B']))

                    proc.stdin.write( "loss_values_data=T\n".encode('latin1'))
                elif k not in transl:
                    continue
                else:
                    proc.stdin.write( "{}=".format(transl[k]).encode('latin1'))
                    if type(mcv[k]) is list:
                        proc.stdin.write( ','.join(map(str,mcv[k])).encode('latin1'))
                    else:
                        proc.stdin.write( "{}".format(mcv[k]).encode('latin1'))
                    proc.stdin.write( "\n".encode('latin1') )

            proc.stdin.write("/\n".encode('latin1'))
            proc.stdin.close()
        except OSError as e:
            logger.error("PATH {}\n CWD {}\n Prog {}\n{}\n".format(
                os.environ['PATH'], directory, writeproc, str(e)))
            e.args += writeproc
            raise e
        except IOError as e:
            logger.error("{} {}\n".format(mcv['name'], str(e)))
            return None
        except KeyError as e:
            logger.error("{} key {} not found\n".format(mcv['name'], str(e)))
            return None
        #for l in proc.stdout:
        #    print l.strip()

        for l in proc.stderr:
            logger.error( l )
        
        return filename
    
    def fitLossCoeffs( self ):
        for m in self.mcv:
            if 'losses' not in self.mcv[m]: continue
            losses = self.mcv[m]['losses']
            # make block matrix
            maxlen=0
            for p in losses['pfe']:
                 if len(p)>maxlen: maxlen=len(p)
            for p in losses['pfe']:
                if len(p)<maxlen: p += [None]*(maxlen-len(p))
            cw,alfa,beta = femagtools.losscoeffs.fit( losses['f'],
                                                      losses['B'],
                                                      losses['pfe'],
                                                      self.mcv[m]['Bo'],
                                                      self.mcv[m]['fo'] )
            losses['cw'] = cw
            losses['alfa'] = alfa
            losses['beta'] = beta
            losses['Bo'] = self.mcv[m]['Bo']
            losses['fo'] = self.mcv[m]['fo']
            
if __name__ == "__main__":
    if len(sys.argv)==2 :
        filename=sys.argv[1]
    else:
        filename=sys.stdin.readline().strip()

    file=open(filename)
    mcvdata=json.load(file)
    file.close()
        
    mcv = MagnetizingCurve(mcvdata)
    mcv.writefile(mcvdata['name'], '.','cat')
    #mcv.recalc()
    #mcv.writefile(mcvdata['name'], '.','cat')
