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
from six import string_types
import femagtools.losscoeffs

# curve types
MAGCRV = 1
DEMCRV = 2
ORIENT_CRV = 3
MAG_AC_CRV = -1
ORIENT_PM_CRV = 5
DEMCRV_BR = 4

logger = logging.getLogger(__name__)

transl = dict(
    cversion='version_mc_curve',
    desc='mc1_title',

    ni='mc1_ni',
    mi='mc1_mi',

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
    curve=[dict(
        hi='mc1_hi',
        bi='mc1_bi',
        bi2='mc1_bi2',
        nuer='mc1_nuer',
        a='mc1_a',
        b='mc1_b',
    )],
    db2='mc1_db2',
    fe_sat_mag='mc1_fe_sat_magnetization'
)

MC1_MIMAX = 50
MC1_NIMAX = 50
M_LOSS_INDUCT = 20
M_LOSS_FREQ = 20

MUE0 = 4e-7*np.pi  # 1.2566371E-06


def approx(db2, curve):
    """return nuer, bi2, a, b approx for curve"""
    nuek0 = (curve['hi'][1] - curve['hi'][0]) / \
            (curve['bi'][1]-curve['bi'][0])
    bk02 = curve['bi'][1]**2
    nuer = [MUE0*nuek0]
    bi2 = [bk02]
    a = []
    b = []

    bk1 = 0.0
    while bk1 <= curve['bi'][-1]:
        bk12 = bk02 + db2
        bk1 = np.sqrt(bk12)
        j = 1
        while j < len(curve['bi']) and bk1 > curve['bi'][j]:
            j += 1
        j -= 1
        bdel = curve['bi'][j] - curve['bi'][j-1]
        c1 = (curve['hi'][j] - curve['hi'][j-1])/bdel
        c2 = curve['hi'][j-1] - c1*curve['bi'][j-1]

        nuek1 = c1 + c2/bk1
        a.append(MUE0*(bk12*nuek0 -
                       bk02*nuek1)/db2)
        b.append(MUE0*(nuek1 - nuek0)/db2)
        nuek0 = nuek1
        bk02 = bk12

        nuer.append(MUE0*nuek1)
        bi2.append(bk12)

    a.append(1.0)
    b.append(MUE0*curve['hi'][-1]-curve['bi'][-1])
    return dict(nuer=nuer, a=a, b=b, bi2=bi2)

    
def findNotNone(l):
    """return lower and upper indexes of not none values in list"""
    for i in range(len(l)):
        if l[i]:
            break
    for j in range(len(l)-1, -1, -1):
        if l[j]:
            break
    return (i, j)


class Mcv:
    def __init__(self):
        # default values from file: mcv.par
        self.ACT_VERSION_MC_CURVE = 0
        self.ORIENTED_VERSION_MC_CURVE = 1
        self.PARAMETER_PM_CURVE = 2

        self.MC1_BASE_FREQUENCY = 50.0
        self.MC1_BASE_INDUCTION = 1.5
        self.MC1_CH_FACTOR = 0.0
        self.MC1_CW_FACTOR = 0.0
        self.MC1_CH_FREQ_FACTOR = 0.0
        self.MC1_CW_FREQ_FACTOR = 0.0
        self.MC1_INDUCTION_FACTOR = 0.0
        self.MC1_FE_SPEZ_WEIGTH = 7.65
        self.MC1_FE_SAT_MAGNETIZATION = 2.15

        self.mc1_base_frequency = self.MC1_BASE_FREQUENCY
        self.mc1_base_induction = self.MC1_BASE_INDUCTION
        self.mc1_ch_factor = self.MC1_CH_FACTOR
        self.mc1_cw_factor = self.MC1_CW_FACTOR
        self.mc1_ch_freq_factor = self.MC1_CH_FREQ_FACTOR
        self.mc1_cw_freq_factor = self.MC1_CW_FREQ_FACTOR
        self.mc1_induction_factor = self.MC1_INDUCTION_FACTOR
        self.mc1_fe_spez_weigth = self.MC1_FE_SPEZ_WEIGTH
        self.mc1_fe_sat_magnetization = self.MC1_FE_SAT_MAGNETIZATION
        
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

        self.curve = []
        self.mc1_mi = [self.MC1_MIMAX-2]*self.MCURVES_MAX
        self.mc1_db2 = [0.]*self.MCURVES_MAX
        self.mc1_angle = [0.]*self.MCURVES_MAX

        self.mc1_ni = [0]*self.MCURVES_MAX

        self.mc1_energy = [[0]*self.MCURVES_MAX]*self.MC1_NIMAX

    def rtrimValueList(self, vlist):
        """cut list at first 0"""
        le = len(vlist)
        for i in range(le-1, -1, -1):
            if vlist[i] != 0.:
                break
        return list(vlist[:i+1])


class Writer(Mcv):
    def __init__(self, data=None):
        Mcv.__init__(self)
        if data:
            self.setData(data)

    def __setattr__(self, name, val):
        try:
            self.__dict__[name] = val
        except:  # file format unknown
            logger.debug("setAttr Exception, name: %s, value: %s", name, val)

    def setData(self, data):
        wtrans = {transl[k]: k
                  for k in transl if not isinstance(transl[k], list)}
        for k in wtrans:
            if wtrans[k] in data:
                self.__setattr__(k, data[wtrans[k]])
        self.curve = [dict(bi=c['bi'], hi=c['hi'])
                      for c in data['curve']]
        try:
            self.mc1_angle = [c['angle'] for c in data['curve']]
        except:
            pass
        try:
            self.losses = data['losses']
        except:
            pass
        return

    def getBlockLength(self, d):
        if isinstance(d, string_types) or isinstance(d, bytes):
            try:
                s = bytes(d).decode('utf-8').encode('latin1')
            except:
                s = d.encode('latin1')
            return len(s)
        elif isinstance(d, int) or isinstance(d, float):
            return 4
        elif isinstance(d, list):
            le = 4 * len(d)
            if len(d) and isinstance(d[0], tuple):
                le *= len(d[0])
            return le
        elif isinstance(d, zip):
            import copy
            dc = copy.deepcopy(d)
            return 4 * sum(len(i) for i in dc)
        return None

    def writeBlock(self, d):
        le = self.getBlockLength(d)
        self.fp.write(struct.pack('i', le))
        if isinstance(d, string_types):
            try:
                s = bytes(d).decode('utf-8').encode('latin1')
            except:
                s = d.encode('latin1')
            self.fp.write(s)
        elif isinstance(d, int):
            self.fp.write(struct.pack('i', d))
        elif isinstance(d, float):
            self.fp.write(struct.pack('f', d))
        elif isinstance(d, list):
            for i in d:
                self.writeData(i)
        elif isinstance(d, zip):  # python3
            for i in d:
                self.writeData(i)
        else:
            pass
        self.fp.write(struct.pack('i', le))

    def writeData(self, d):
        if isinstance(d, string_types):
            self.fp.write(bytes(d).decode('utf-8').encode('latin1'))
        elif isinstance(d, int):
            self.fp.write(struct.pack('i', d))
        elif isinstance(d, tuple):
            for i in d:
                self.writeData(i)
        else:
            # must be float?
            self.fp.write(struct.pack('f', d))

    def _prepare(self):
        """prepare output format (internal use only)"""
        self.mc1_curves = len(self.curve)
        self.mc1_ni = [min(len(c['hi']),
                           len(c['bi']))
                       for c in self.curve if 'hi' in c]
        self.mc1_db2 = [(c['bi'][-1]**2 - c['bi'][0]**2)/n
                        for c, n in zip(self.curve, self.mc1_mi)]
        for db2, c in zip(self.mc1_db2, self.curve):
            c.update(approx(db2, c))
        self.mc1_mi = [len(c['a'])
                       for c in self.curve]
        
    def writeBinaryFile(self):
        self._prepare()
        # write line, version_mc_curve
        self.writeBlock(self.version_mc_curve)

        # write line, text '    *** File with magnetic curve ***    '
        self.writeBlock('    *** File with magnetic curve ***    ')
                    
        # write line, mc1_title
        self.writeBlock(self.mc1_title.ljust(40))
        # write line, mc1_ni(1),mc1_mi(1),mc1_type,mc1_recalc,mc1_db2(1)
        self.writeBlock([int(self.mc1_ni[0]),
                         int(self.mc1_mi[0]),
                         int(self.mc1_type),
                         int(self.mc1_recalc),
                         self.mc1_db2[0]])

        # write line, mc1_remz, mc1_bsat, mc1_bref, mc1_fillfac
        if self.version_mc_curve == self.ACT_VERSION_MC_CURVE:
            self.writeBlock([self.mc1_remz, self.mc1_bsat,
                             self.mc1_bref, self.mc1_fillfac])
        if self.mc1_type == DEMCRV_BR:
            self.mc1_remz = self.mc1_angle[self.mc1_curves]
        if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
           self.version_mc_curve == self.PARAMETER_PM_CURVE:
            self.writeBlock([self.mc1_remz, self.mc1_bsat,
                             self.mc1_bref, self.mc1_fillfac,
                             self.mc1_curves])

        if self.mc1_type == DEMCRV_BR:
            self.mc1_angle[self.mc1_curves] = self.mc1_remz

        # data
        for K in range(0, self.mc1_curves):

            # hi, bi
            lb = self.curve[K].get('bi', [])
            lh = self.curve[K].get('hi', [])
            self.writeBlock(zip(*[
                [float(lb[j]) if j < len(lb) else 0.
                 for j in range(self.MC1_NIMAX)],
                [float(lh[j]) if j < len(lh) else 0.
                 for j in range(self.MC1_NIMAX)]]))

            # bi2, nuer
            lb = self.curve[K]['bi2']
            ln = self.curve[K]['nuer']
            self.writeBlock(zip(*[
                [float(lb[j]) if j < len(lb) else 0.
                 for j in range(self.MC1_NIMAX)],
                [float(ln[j]) if j < len(ln) else 0.
                 for j in range(self.MC1_NIMAX)]]))

            # a, b, c, d
            la = self.curve[K].get('a', [0.]*self.MC1_NIMAX)
            lb = self.curve[K].get('b', [0.]*self.MC1_NIMAX)
            self.writeBlock(zip(*[
                [float(la[j]) if j < len(la) else 0.
                 for j in range(self.MC1_NIMAX)],
                [float(lb[j]) if j < len(lb) else 0.
                 for j in range(self.MC1_NIMAX)],
                [0.]*50,
                [0.]*50
            ]))

            #
            if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
               self.version_mc_curve == self.PARAMETER_PM_CURVE:
                self.writeBlock([self.mc1_angle[K], self.mc1_db2[K]])

        self.writeBlock([float(self.mc1_base_frequency),
                         float(self.mc1_base_induction),
                         float(self.mc1_ch_factor),
                         float(self.mc1_cw_factor),
                         float(self.mc1_ch_freq_factor),
                         float(self.mc1_cw_freq_factor),
                         float(self.mc1_induction_factor),
                         float(self.mc1_fe_spez_weigth),
                         float(self.mc1_fe_sat_magnetization)])

        try:
            nfreq = len([1 for x in self.losses['f'] if x > 0])
            nind = len(self.losses['B'])
            self.writeBlock([nfreq, nind])
            self.writeBlock(self.losses['B'] +
                            [0.0]*(M_LOSS_INDUCT - len(self.losses['B'])))
            mloss = M_LOSS_FREQ
            cw = self.mc1_cw_factor
            alpha = self.mc1_cw_freq_factor
            ch = self.mc1_ch_factor
            beta = self.mc1_ch_freq_factor
            gamma = self.mc1_induction_factor
            jordan = lambda fr, Br: (
                (cw*fr**alpha + ch*fr**beta) * Br**gamma)
            i = 1
            for f, p in zip(self.losses['f'], zip(*self.losses['pfe'])):
                if f:
                    pl = [px if px else jordan(f/self.losses['fo'],
                                               b/self.losses['Bo'])
                          for px, b in zip(p, self.losses['B'])]

                    self.writeBlock(pl +
                                    [0.0]*(M_LOSS_INDUCT - len(p)))
                    self.writeBlock([f])
                    i += 1
                    mloss -= 1
            for m in range(mloss):
                self.writeBlock([0.0]*M_LOSS_INDUCT)
                self.writeBlock([0.0])
                i += 1

            self.writeBlock([self.losses['cw'], self.losses['cw_freq'],
                             self.losses['b_coeff'], self.losses['Bo'],
                             self.losses['fo']])
            self.writeBlock([1])
            logger.info('Losses n freq %d n ind %d', nfreq, nind)
        except:
            pass
            
    def writeMcv(self, filename):
        # windows needs this strip to remove '\r'
        filename = filename.strip()
        self.name = os.path.splitext(filename)[0]

        if filename.upper().endswith('.MCV') or \
           filename.upper().endswith('.MC'):
            binary = True
            self.fp = open(filename, "wb")
        else:
            binary = False
            self.fp = open(filename, "wb")
        logger.info("Write File %s, binary format %d", filename, binary)

        self.writeBinaryFile()
        self.fp.close()


class Reader(Mcv):
    def __init__(self):
        Mcv.__init__(self)

    def readBlock(self, d, length=0):
        res = []
        try:
            le = self.getInteger()
            if d == string_types:
                res = self.getString(length)
            elif d == int:
                res = self.getInteger()
            elif d == float:
                res = self.getReal()
            elif isinstance(d, list):
                res = [self.readData(i) for i in d]
            else:
                pass

            le2 = self.getInteger()
        except:  # file format unknown
            le2 = 0
#        logger.debug("readBlock Len: %d == %d, data: %s", le, le2, res)
        return res

    def readData(self, d, length=1):
        if d == string_types:
            return self.getString(length)
        elif d == int:
            return self.getInteger()
        else:
            # must be float?
            return self.getReal()

    def getString(self, length=1):
        block = self.fp.read(length)
        st = []
        for i in range(0, len(block)):
            (s,) = struct.unpack('c', block[i:i+1])
            st.append(s.decode('latin1'))
        return ''.join(st)

    def getInteger(self, length=4):
        block = self.fp.read(length)
        (integer,) = struct.unpack('i', block[0:length])
        return integer

    def getReal(self):
        block = self.fp.read(4)
        if len(block) == 4:
            (real,) = struct.unpack('f', block[0:4])
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

        self.name = os.path.splitext(os.path.basename(filename))[0]
        # read curve version (INTEGER)
        if binary:
            self.version_mc_curve = self.readBlock(int)
        else:
            self.version_mc_curve = int(self.fp.readline().strip())
        logger.info("MC Version %s", self.version_mc_curve)

        # read dummy text and title 2x (CHARACTER*40)
        if binary:
            # info text '*** File with magnetic curve ***'
            self.readBlock(string_types, 40)
            self.mc1_title = self.readBlock(string_types, 40)  # read title
        else:
            self.fp.readline().strip()
            self.mc1_title = self.fp.readline().strip()

        # read line 4
        if binary:
            (self.mc1_ni[0],
             self.mc1_mi[0],
             self.mc1_type,
             self.mc1_recalc,
             self.mc1_db2[0]) = self.readBlock([int, int, int, int, float])
        else:
            line = self.fp.readline().split()
            self.mc1_ni[0]  = int (line[0])   # mc1_ni (INTEGER)
            self.mc1_mi[0]  = int (line[1])   # mc1_mi (INTEGER)
            self.mc1_type   = int (line[2])   # mc1_type (INTEGER)
            self.mc1_recalc = int (line[3])   # mc1_recalc (INTEGER)
            self.mc1_db2[0] = float (line[4]) # mc1_db2 (REAL)

        # read line 5
        if binary:
            l_format = [float, float, float, float]
            if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
               self.version_mc_curve == self.PARAMETER_PM_CURVE:
                l_format.append(int)
            t = self.readBlock(l_format)
            self.mc1_remz = t[0]
            self.mc1_bsat = t[1]
            self.mc1_bref = t[2]
            self.mc1_fillfac = t[3]
            if len(t) > 4:
                self.mc1_curves = t[4]
        else:
            line = self.fp.readline()
            (self.mc1_remz, self.mc1_bsat, self.mc1_bref, self.mc1_fillfac) = \
                            map(float, line.split())

            if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
               self.version_mc_curve == self.PARAMETER_PM_CURVE:
                self.mc1_curves = int(line[4])

        if self.mc1_type == DEMCRV_BR:
            self.mc1_angle[self.mc1_curves] = self.mc1_remz

        if not binary:
            # read rest of file and convert all to float values
            values = map(float, ' '.join(self.fp.readlines()).split())

        self.curve = []
        for K in range(0, self.mc1_curves):
            if binary:
                # bi, hi
                res = self.readBlock([float]*2*self.MC1_MIMAX)
                mc_bi = res[::2]
                mc_hi = res[1::2]

                # bi2, nuer
                res = self.readBlock([float]*2*self.MC1_MIMAX)
                mc_bi2 = res[::2]
                mc_nuer = res[1::2]

                # a, b, c, d
                res = self.readBlock([float]*4*self.MC1_MIMAX)
                mc_a = res[::4]
                mc_b = res[1::4]
                mc_c = res[2::4]
                mc_d = res[3::4]
            else:
                [mc_bi, mc_hi] = zip(*[values[2*I:2*I+2]
                                       for I in range(self.MC1_NIMAX)])
                idxOffset = 2*self.MC1_NIMAX
                [mc_bi2, mc_nuer] = zip(*[values[idxOffset+2*I:idxOffset+2*I+2]
                                          for I in range(self.MC1_NIMAX)])
                idxOffset += 2*self.MC1_NIMAX
                [mc_a, mc_b, mc_c, mc_d] = zip(*[
                    values[idxOffset+4*I:idxOffset+4*I+4]
                    for I in range(self.MC1_NIMAX)])
                idxOffset += 4*self.MC1_NIMAX

            if self.version_mc_curve == self.ORIENTED_VERSION_MC_CURVE or \
               self.version_mc_curve == self.PARAMETER_PM_CURVE:
                if binary:
                    (self.mc1_angle[K],
                     self.mc1_db2[K]) = self.readBlock([float]*2)
                else:
                    (self.mc1_angle[K],
                     self.mc1_db2[K]) = (values[idxOffset:idxOffset+2])
                    idxOffset += 2

            for I in range(self.MC1_NIMAX):
                if mc_bi[I] != 0.0 or mc_hi[I] != 0.0:
                    self.mc1_ni[K] = I+1
            for I in range(self.MC1_MIMAX):
                if mc_a[I] != 0.0 or mc_b[I] != 0.0:
                    self.mc1_mi[K] = I+1
            # assign data
            self.curve.append(dict(
                hi=self.rtrimValueList(
                    [mc_hi[I] for I in range(self.MC1_NIMAX)]),
                bi=self.rtrimValueList(
                    [mc_bi[I] for I in range(self.MC1_NIMAX)]),
                bi2=self.rtrimValueList(
                    [mc_bi2[I] for I in range(self.MC1_MIMAX)]),
                nuer=self.rtrimValueList(
                    [mc_nuer[I] for I in range(self.MC1_MIMAX)]),
                a=self.rtrimValueList(
                    [mc_a[I] for I in range(self.MC1_MIMAX)]),
                b=self.rtrimValueList(
                    [mc_b[I] for I in range(self.MC1_MIMAX)])))

        # set dummy defaults
        vals = [self.MC1_BASE_FREQUENCY,
                self.MC1_BASE_INDUCTION,
                self.MC1_CH_FACTOR,
                self.MC1_CW_FACTOR,
                self.MC1_CH_FREQ_FACTOR,
                self.MC1_CW_FREQ_FACTOR,
                self.MC1_INDUCTION_FACTOR,
                self.MC1_FE_SPEZ_WEIGTH,
                self.MC1_FE_SAT_MAGNETIZATION]

        if binary:
            res = self.readBlock([float]*9)
            for i in range(len(res)):
                if math.isnan(res[i]):
                    break
                vals[i] = res[i]
        else:
            iLen = min(int(s) for s in [9, (len(values)-idxOffset)])
            for I in range(iLen):
                vals[I] = values[idxOffset+I]
            idxOffset += iLen
        logger.debug("Mcv last Props: %s", vals)

        self.fo = vals[0]
        self.Bo = vals[1]
        self.ch = vals[2]
        self.ch_freq = vals[4]
        self.cw = vals[3]
        self.cw_freq = vals[5]
        self.b_coeff = vals[6]
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
                
            self.mc1_energy[K] = [energy[I] for I in range(len(energy))]

    def get_results(self):
        result = {
            'name': self.name,
            'desc': self.mc1_title,
            'cversion': self.version_mc_curve,
            #'ni': [n for n in self.mc1_ni if n],
            #'mi': [m for m in self.mc1_mi if m],
            'ctype': self.mc1_type,
            'recalc': self.mc1_recalc,
            #'db2': self.mc1_db2,
            'remz': self.mc1_remz,
            'bsat': self.mc1_bsat,
            'bref': self.mc1_bref,
            'fillfac': self.mc1_fillfac,
            #'curves': self.mc1_curves,
            #'curve': self.curve,
            #'energy': self.mc1_energy,
            'fo': self.fo,
            'Bo': self.Bo,
            'ch': self.ch,
            'ch_freq': self.ch_freq,
            'cw': self.cw,
            'cw_freq': self.cw_freq,
            'b_coeff': self.b_coeff,
            'rho': self.rho,
            'fe_sat_mag': self.fe_sat_mag,
            'curve': [{
                'bi': c.get('bi'),
                'hi': c.get('hi'),
            } for c in self.curve]
        }
        try:
            result['losses'] = self.losses
        except:
            pass
        if (self.ORIENTED_VERSION_MC_CURVE or
            self.PARAMETER_PM_CURVE):
            for i in range(len(self.curve)):
                result['curve'][i]['angle'] = self.mc1_angle[i]
            
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
                return
            if 'name' in mcvpar:
                self.mcv[mcvpar['name']] = mcvpar
                return
            self.mcv['0'] = mcvpar

        # Do not use unicode in PYTHON 3 as all strings are sequences of Unicode
        elif isinstance(mcvpar, string_types) or isinstance(mcvpar, unicode):
            self.mcdirectory = os.path.abspath(mcvpar)
            logger.info("MC Dir %s", self.mcdirectory)
        else:
            raise Exception("unsupported parameter type "+str(type(mcvpar)))

    def find(self, id):
        """find mcv by id or name"""
        try:
            return self.mcv[id]['name']
        except ValueError:
            pass  # not found
        except KeyError:
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
        try:
            for k in self.mcv.keys():
                if self.mcv[k]['name'] == name:
                    return self.mcv[k]
        except KeyError:
            pass
        # not found
        if len(self.mcv) == 1 and '0' in self.mcv:
            self.mcv['0']['name'] = name
            return self.mcv['0']
        return None

    def recalc(self):
        for m in self.mcv:
            curve = self.mcv[m]['curve'][0]
            mi = MC1_MIMAX-2
            dh = curve['hi'][-1]-curve['hi'][-2]
            db = curve['bi'][-1]-curve['bi'][-2]
            dmue_d = db/dh
            dmue = curve['bi'][-1]/curve['hi'][-1]
            db = 3e-2*curve['bi'][-1]
            n3 = 1.5

            curve['muer'] = [b/MUE0/h
                             for b, h in zip(curve['bi'],
                                             curve['hi'])]
            
            # extend bh-curve into saturation
            while dmue_d > 1.01*MUE0 and dmue > 1.5*MUE0:
                dmue_d = MUE0 + (dmue_d-MUE0)/np.sqrt(n3)
                curve['bi'].append(curve['bi'][-1]+db)
                curve['hi'].append(curve['hi'][-1]+db/dmue_d)
                curve['muer'].append(curve['bi'][-1]/MUE0/curve['hi'][-1])
                n3 += 0.2
                dmue = curve['bi'][-1]/curve['hi'][-1]

            self.mcv[m]['db2'] = (curve['bi'][-1]**2 -
                                  curve['bi'][0]**2)/(mi-1)
            nuek0 = (curve['hi'][1] - curve['hi'][0]) / \
                    (curve['bi'][1]-curve['bi'][0])
            for j1 in range(len(curve['bi'])):
                    bk02 = curve['bi'][j1]**2
                    if bk02 > 0:
                        break
            curve['nuer'] = [MUE0*nuek0]
            curve['bi2'] = [bk02]
            curve['a'] = []
            curve['b'] = []

            bk1 = 0.0
            while bk1 <= curve['bi'][-1]:
                    bk12 = bk02 + self.mcv[m]['db2']
                    bk1 = np.sqrt(bk12)
                    j = 2
                    while j < len(curve['bi']) and bk1 <= curve['bi'][j]:
                        j += 1
                    j -= 1
                    bdel = curve['bi'][j] - curve['bi'][j1]
                    c1 = (curve['hi'][j] - curve['hi'][j1])/bdel
                    c2 = curve['hi'][j1] - c1*curve['bi'][j1]

                    nuek1 = c1 + c2/bk1

                    curve['a'].append(MUE0*(bk12*nuek0 -
                                            bk02*nuek1)/self.mcv[m]['db2'])
                    curve['b'].append(MUE0*(nuek1 - nuek0)/self.mcv[m]['db2'])
                    nuek0 = nuek1
                    bk02 = bk12

                    curve['nuer'].append(MUE0*nuek1)
                    curve['bi2'].append(bk12)

            curve['a'].append(1.0)
            curve['b'].append(MUE0*curve['hi'][-1]-curve['bi'][-1])
            
    def writefile(self, name, directory='.', writeproc=''):
        """find magnetic curve by name or id and write binary file
        returns filename if found else None"""
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
            except shutil.SameFileError:
                return filename
            except:
                logger.error("MCV %s not found", str(filename))
            return None

        filename = ''.join((mcv['name'], ext))
        if writeproc:
            self._writefile(mcv, directory, filename, writeproc)
        else:
            writer = Writer(mcv)
            writer.writeMcv(os.path.join(directory, filename))
        return filename

    def _writefile(self, mcv, directory, filename, writeproc):
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
                    pfeT = []
#                    flen=10
#                    blen=20
                    cw = mcv[k]['cw']
                    alfa = mcv[k]['alfa']
                    beta = mcv[k]['beta']
                    fo  = mcv[k]['fo']
                    Bo = mcv[k]['Bo']

                    pfe = []
                    lower = 0
                    # must replace all None by approx losses
                    for i in range(len(mcv['losses']['f'])):
                        f = mcv['losses']['f'][i]
                        if f > 0:
                            pfei = [p[i] if i < len(p) else None
                                    for p in mcv['losses']['pfe']]
                            m, n = findNotNone(pfei)
                            if m > lower:
                                lower = m
                            if m <= n:
                                y = [np.log10(p) for p in pfei[m:n+1]]
                                x = [np.log10(b/Bo)
                                     for b in mcv['losses']['B'][m:n+1]]
                                A = np.vstack([x, np.ones(len(x))]).T
                                beta, cw = np.linalg.lstsq(A, y)[0]
                                for j in range(n+1, len(pfei)):
                                    pfei[j] = 10**cw*(
                                        mcv['losses']['B'][j]/Bo)**beta

                                pfe.append(pfei)
                    pfe += [[0] * M_LOSS_INDUCT] * (
                        M_LOSS_FREQ - len(mcv['losses']['f']))

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
            
    def fitLossCoeffs(self):
        for m in self.mcv:
            if 'losses' not in self.mcv[m]:
                continue
            losses = self.mcv[m]['losses']
            cw, alfa, beta = femagtools.losscoeffs.fitsteinmetz(
                losses['f'],
                losses['B'],
                losses['pfe'],
                self.mcv[m]['Bo'],
                self.mcv[m]['fo'])
            losses['cw'] = cw
            losses['alfa'] = alfa
            losses['beta'] = beta
            losses['Bo'] = self.mcv[m]['Bo']
            losses['fo'] = self.mcv[m]['fo']
            
if __name__ == "__main__":
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    mcv = Reader()
    mcv.readMcv(filename)
    json.dump(mcv.get_results(), sys.stdout)
    #print(mcv.get_results()
    #mcv.recalc()
    #mcv.writefile(mcvdata['name'], '.','cat')
