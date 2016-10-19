# -*- coding: utf-8 -*-
"""
    femagtools.mcv
    ~~~~~~~~~~~~~~

    Reading, Creating and managing MCV/MC files

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import json
import subprocess
import sys
import logging
import os.path
import struct
import math
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


class Reader:
    def __init__(self):
        # default values from file: mcv.par
        self.ACT_VERSION_MC_CURVE = 0
        self.ORIENTED_VERSION_MC_CURVE = 1
        self.PARAMETER_PM_CURVE = 2
        self.DEMCRV_BR = 4

        self.MC1_BASE_FREQUENCY = 50.0
        self.MC1_BASE_INDUCTION = 1.5
        self.MC1_CH_FACTOR = 0.0
        self.MC1_CW_FACTOR = 0.0
        self.MC1_CH_FREQ_FACTOR = 0.0
        self.MC1_CW_FREQ_FACTOR = 0.0
        self.MC1_INDUCTION_FACTOR = 0.0
        self.MC1_FE_SPEZ_WEIGTH = 7.65
        self.MC1_FE_SAT_MAGNETIZATION = 2.15

        self.mc1_title = ''
        self.version_mc_curve = self.ACT_VERSION_MC_CURVE
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

        self.curve=[]
        self.mc1_mi    = [ 0. for I in range(self.MCURVES_MAX)  ]
        self.mc1_db2   = [ 0. for I in range(self.MCURVES_MAX)  ]
        self.mc1_angle = [ 0. for I in range(self.MCURVES_MAX)  ]

        self.mc1_ni = [ float('nan') for I in range(self.MCURVES_MAX)  ]

        self.mc1_bi=[ [ float('nan') for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_NIMAX) ]
        self.mc1_hi=[ [ 0. for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_NIMAX) ]

        self.mc1_bi2=[ [ 0. for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_MIMAX) ]
        self.mc1_nuer=[ [ 0. for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_MIMAX) ]
        self.mc1_a=[ [ 0. for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_MIMAX) ]
        self.mc1_b=[ [ 0. for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_MIMAX) ]
        self.mc1_energy=[ [ 0 for I in range(self.MCURVES_MAX)  ] for K in range(self.MC1_MIMAX) ]

    def getString(self, length=1):
        block = self.fp.read(length)
        st = []
        for i in range(0, len(block)):
            (s,) = struct.unpack('c', block[i:i+1])
            if ord(s) != 4 and ord(s) != 0  and ord(s) != 12  and ord(s) != 16 :
                st.append(s)
#        print "getString: [", "", st, "][", "", st.strip(), "] REAL[",real,"] len(", len(st), ') lenStrip(', len(st.strip()), ')\n'
        return b''.join(st).decode('latin1').strip().encode('utf-8')
        
    def getInteger(self):
        block = self.fp.read(4)
        (integer,) = struct.unpack('i', block[0:4])
#        print "Integer: ", integer
        return integer

    def rtrimValueList(self, vlist):
        le = len(vlist)
        for i in range(le-1, -1, -1):
            if vlist[i] == 0.:
                vlist = vlist[:-1]
            else:
                break
        return vlist
  
    def getReal(self):
        block = self.fp.read(4)
        if len(block) == 4 and ord(block[0:1]) in [12, 16] and \
               ord(block[1:2]) == 0 and ord(block[2:3]) == 0 and \
               ord(block[3:4]) == 0:
#            print "Read Next "
            return self.getReal()
#        for i in range(0, len(block)):
#            (s,) = struct.unpack('c', block[i:i+1])
#            print "RE C:", ord(s)
        if len(block) == 4:
            (real,) = struct.unpack('f', block[0:4])
#            print "Real: ", real
            return real
        else:
            return float('nan')
 
    def readMcv(self, filename):
        # intens bug : windows needs this strip to remove '\r'
        filename = filename.strip()

        if filename.upper().endswith('.MCV') or \
           filename.upper().endswith('.MC'):
            binary = True
            self.fp = open(filename, "rb")
        else:
            binary = False
            self.fp = open(filename, "r")

        self.name = os.path.splitext(filename)[0]
        # read curve version (INTEGER)
        if binary:
            str = self.getString(4)  # dummy 4 '\4\0\0\0'
            self.version_mc_curve = self.getInteger()
            str = self.getString(8)  # dummy 8 '\4\0\0\0' + '(' + '\0\0\0'
        else:
            self.version_mc_curve = int(self.fp.readline().strip())
    
        # read dummy text and title 2x (CHARACTER*40)
        if binary:
            str = self.getString(40)  # info text '*** File with magnetic curve ***'
            str = self.getString(8)  # dummy 8 '(('
            self.mc1_title = self.getString(40)  # read title
        else:
            self.fp.readline().strip()
            self.mc1_title = self.fp.readline().strip()
            
        #self.mc1_title = self.mc1_title.encode('utf-8')

        #read line 4
        if binary:
            str = self.getString(8) # dummy 8 '('+'\0\0\0\20\0\0\0'
            self.mc1_ni[0]  = self.getInteger()  # mc1_ni (INTEGER)
            self.mc1_mi[0]  = self.getInteger()  # mc1_mi (INTEGER)
            self.mc1_type   = self.getInteger()  # mc1_type (INTEGER)
            self.mc1_recalc = self.getInteger()  # mc1_recalc (INTEGER)
            self.mc1_db2[0] = self.getReal()     # mc1_db2 (REAL)
        else:
            line = self.fp.readline().split()
            self.mc1_ni[0]  = int (line[0])   # mc1_ni (INTEGER)
            self.mc1_mi[0]  = int (line[1])   # mc1_mi (INTEGER)
            self.mc1_type   = int (line[2])   # mc1_type (INTEGER)
            self.mc1_recalc = int (line[3])   # mc1_recalc (INTEGER)
            self.mc1_db2[0] = float (line[4]) # mc1_db2 (REAL)

        #read line 5
        if binary:
            self.getString(8) # dummy 8
            self.mc1_remz = self.getReal()
            self.mc1_bsat = self.getReal()
            self.mc1_bref = self.getReal()
            self.mc1_fillfac = self.getReal()
            if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
                   self.version_mc_curve == self.PARAMETER_PM_CURVE:
                self.mc1_curves = self.getInteger()

        else:
            line = self.fp.readline()
            (self.mc1_remz, self.mc1_bsat, self.mc1_bref, self.mc1_fillfac) = \
                            map( float, line.split())
            
        if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
               self.version_mc_curve == self.PARAMETER_PM_CURVE:
            if binary:
                self.mc1_curves = self.getInteger()
            else:
                self.mc1_curves = int(line[4])

        if self.mc1_type == self.DEMCRV_BR :
            self.mc1_angle[mc1_curves] = self.mc1_remz

        if binary:
            self.getString(8) # dummy 8
        else:
            # read rest of file and convert all to float values
            values = map( float, ' '.join(self.fp.readlines()).split())

        self.curve=[{} for i in range(self.mc1_curves)]
        for K in range(0, self.mc1_curves):
            if binary:
                [mc_bi,mc_hi]=zip(*[(self.getReal(), self.getReal())
                                    for I in range(self.MC1_NIMAX) ])
#                self.getString(8) # dummy 8
                self.getString(8)
                [mc_bi2, mc_nuer]=zip(*[(self.getReal(), self.getReal())
                                       for I in range(self.MC1_MIMAX) ])
                self.getString(8) # dummy 8
                [mc_a,mc_b, mc_c, mc_d]=zip(*[(self.getReal(), self.getReal(),
                                               self.getReal(), self.getReal())
                                              for I in range(self.MC1_MIMAX) ])
                self.getString(8) # dummy 8
            else:
                [mc_bi,mc_hi]= zip(*[values[2*I:2*I+2] for I in range(self.MC1_NIMAX)])
                idxOffset = 2*self.MC1_NIMAX
                [mc_bi2,mc_nuer]= zip(*[values[idxOffset+2*I:idxOffset+2*I+2]
                                        for I in range(self.MC1_NIMAX)])
                idxOffset += 2*self.MC1_NIMAX
                [mc_a,mc_b, mc_c, mc_d]= zip(*[values[idxOffset+4*I:idxOffset+4*I+4]
                                               for I in range(self.MC1_NIMAX)])
                idxOffset += 4*self.MC1_NIMAX

            if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
                   self.version_mc_curve == self.PARAMETER_PM_CURVE :
                if binary:
                    (self.mc1_angle[K], self.mc1_db2[K]) = (self.getReal(), self.getReal())
                else:
                    (self.mc1_angle[K], self.mc1_db2[K])=  \
                                        (values[idxOffset:idxOffset+2])
                    idxOffset += 2


            self.curve[K]['hi'] = self.rtrimValueList(
                [ mc_hi[I] for I in range(self.MC1_NIMAX)  ] )
            self.curve[K]['bi'] = self.rtrimValueList(
                [ mc_bi[I] for I in range(self.MC1_NIMAX)  ] )
            
            for I in range(self.MC1_NIMAX):
                if mc_bi[I] != 0.0 or mc_hi[I] != 0.0 :
                    self.mc1_ni[K] = I+1

            self.curve[K]['bi2']  = self.rtrimValueList(
                [ mc_bi2[I] for I in range(self.MC1_MIMAX)  ] )

            self.curve[K]['nuer'] = self.rtrimValueList(
                [ mc_nuer[I] for I in range(self.MC1_MIMAX)  ] )
            self.curve[K]['a'] = self.rtrimValueList(
                [ mc_a[I] for I in range(self.MC1_MIMAX)  ] )
            self.curve[K]['b'] = self.rtrimValueList(
                [ mc_b[I] for I in range(self.MC1_MIMAX)  ] )
            
            for I in range(self.MC1_MIMAX):
                if mc_a[I] !=  0.0 or mc_b[I] != 0.0 : 
                    self.mc1_mi[K] = I+1

        # set dummy defaults
        vals = [ self.MC1_BASE_FREQUENCY,
                 self.MC1_BASE_INDUCTION,
                 self.MC1_CH_FACTOR,
                 self.MC1_CW_FACTOR,
                 self.MC1_CH_FREQ_FACTOR,
                 self.MC1_CW_FREQ_FACTOR,
                 self.MC1_INDUCTION_FACTOR,
                 self.MC1_FE_SPEZ_WEIGTH,
                 self.MC1_FE_SAT_MAGNETIZATION]

        #print vals

        if binary:
            self.getString(8) # dummy 8
            for I in range(9):
                f = self.getReal()
                if math.isnan(f):
                    break
#                print "set dummy I: ", I, "  value: ", f
                vals[I] = f
        else:
            iLen = min(int(s) for s in [9, (len(values)-idxOffset)])
            for I in range(iLen):
                vals[I] = values[idxOffset+I]
            idxOffset += iLen
     
        #print vals

        self.fo = vals[0]
        self.Bo = vals[1]
        self.ch = vals[2]
        self.ch_freq = vals[4]
        self.cw = vals[3]
        self.cw_freq = vals[5]
        self.b_coeff =  vals[6]
        self.rho = vals[7]
        self.fe_sat_mag = vals[8]

        if self.MC1_INDUCTION_FACTOR > 2.0:
            self.MC1_INDUCTION_FACTOR = 2.0
        for K in range(0, self.mc1_curves):
            energy = []
            if binary:
                for I in range(9):
                    f = self.getReal()
                    if math.isnan(f):
                        break
                    energy.append(f)
            else:
                iLen = min([int(s)
                            for s in [self.MC1_MIMAX,
                                      (len(values)-idxOffset)]])
                for I in range(iLen):
                    energy.append(values[idxOffset+I])
                idxOffset += iLen
                
            self.mc1_energy[K] = [ energy[I] for I in range(len(energy))  ]

    def get_results(self):
#        print "version : ", self.version_mc_curve,
#        print "title : ", self.mc1_title
#        print "hi : ", self.mc1_hi
#        print "bi : ", self.mc1_bi
        result={
            'name':self.name,
            'desc': self.mc1_title.decode('utf-8'),
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
            'fe_sat_mag':self.fe_sat_mag }
      
#        print "Title: [",  self.mc1_title, "] CurveVersion :", self.version_mc_curve
#        print "LINE 0: ", self.mc1_ni[0], ' ',  self.mc1_mi[0], ' ',  self.mc1_type, ' ',  self.mc1_recalc, ' ',  self.mc1_db2[0]
#        print "LINE 1: ", self.mc1_remz, ' ',  self.mc1_bsat, ' ',  self.mc1_bref, ' ',  self.mc1_fillfac, ' ',  self.mc1_curves

        return result


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
            
    def writefile(self, name, directory='.', writeproc='mcvwriter'):
        """find magnetic curve by name or id and write binary file
        returns filename if found else None"""
        if not id:
            return None
        ext = '.MC' if sys.platform == 'win32' else '.MCV'
        mcv = self.find_by_name(name)
        if not mcv:
            filename = ''.join((name, ext))
            try:
                import shutil
                logger.info("Copy file %s", filename)
                shutil.copy(os.path.join(self.mcdirectory,
                                         filename), directory)
                return filename
            except:
                logger.error("MCV %s not found", str(filename))
            return None

        filename = ''.join((mcv['name'], ext))
        try:
            logger.info("create %s", str(filename))
            proc = subprocess.Popen([writeproc],
                                    cwd=directory,
                                    stdin=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    bufsize=1)

            inparams = ["&INPARAM"]
            for k in mcv.keys():
                if k == 'name':
                    inparams.append("filename='{}'".format(filename))
                elif k == 'desc':
                    inparams.append("mc1_title='{}'".format(mcv['desc']))
                elif k == 'curve':
                    inparams.append('mc1_curves={}'.format(len(mcv[k])))
                    for c in mcv['curve']:
                        for n in c.keys():
                            if n not in transl:
                                continue
                            inparams.append("{}=".format(transl[n]))
                            if type(c[n]) is list:
                                inparams.append(
                                    ','.join(map(str, c[n])))
                            else:
                                inparams.append(str(c[n]))

                    if 'energy' in mcv['curve'][0] and len(mcv['curve'][0]['energy']) > 0:
                        inparams.append("mc1_energy=")
                        inparams.append(','.join(map(str, mcv['curve'][0]['energy'])))

                elif k == 'losses':
                    # find start index
                    bstart = 0
                    for a in mcv[k]['pfe']:
                        for i in range(len(a)):
                            if a[i]:
                                if i > bstart:
                                    bstart = i
                                break
                    inparams.append('N_freq={}'.format(len(mcv[k]['f'])))
                    inparams.append('N_J_ind={}'.format(len(mcv[k]['B'])-bstart))
                    inparams.append("frequency=")
                    inparams.append(",".join(map(str, mcv[k]['f'])))
                    inparams.append("induction=")
                    inparams.append(",".join(map(str, mcv[k]['B'][bstart:])))
                    inparams.append("losses=")
                    pfeT= []
#                    flen=10
#                    blen=20
                    cw = mcv[k]['cw']
                    alfa = mcv[k]['alfa']
                    beta = mcv[k]['beta']
                    fo  =mcv[k]['fo']
                    Bo = mcv[k]['Bo']

                    pfe = []
                    lower = 0
                    # must replace all None by approx losses
                    for i in range(len(mcv['losses']['f'])):
                        f = mcv['losses']['f'][i]
                        if f > 0:
                            pfei = [p[i] if i<len(p) else None for p in mcv['losses']['pfe']]
                            m, n = findNotNone(pfei)
                            if m > lower: lower=m
                            if m <= n:
                                y = [ np.log10(p) for p in pfei[m:n+1] ]
                                x = [ np.log10(b/Bo) for b in mcv['losses']['B'][m:n+1]]
                                A = np.vstack([x, np.ones(len(x))]).T
                                beta, cw = np.linalg.lstsq(A, y)[0]
                                for j in range(n+1,len(pfei)):
                                    pfei[j] = 10**cw*(mcv['losses']['B'][j]/Bo)**beta

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
                        
                    inparams.append(','.join(map(str, pfeT)))
                    inparams.append("cw_m={}".format(mcv[k]['cw']))
                    inparams.append("alfa_m={}".format(mcv[k]['alfa']))
                    inparams.append("beta_m={}".format(mcv[k]['beta']))
                    inparams.append("base_freq={}".format(mcv[k]['fo']))
                    inparams.append("base_induct={}".format(mcv[k]['Bo']))
                    logger.info("has losses N_freq %d, N_J_ind %d",
                                len(mcv[k]['f']), len(mcv[k]['B']))

                    inparams.append("loss_values_data=T")
                elif k not in transl:
                    continue
                else:
                    inparams.append("{}=".format(transl[k]))
                    if type(mcv[k]) is list:
                        inparams.append(','.join(map(str, mcv[k])))
                    else:
                        inparams.append("{}".format(mcv[k]))

            inparams.append("/\n")
            proc.stdin.write(bytes('\n'.join(inparams).encode('latin-1')))
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

        for l in proc.stderr:
            logger.error(l)
        
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
