#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
import sys
import struct
import math
import json
import codecs
import os

class McvReader:
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
 
    def readMcv( self, filename):
        # intens bug : windows needs this strip to remove '\r'
        filename = filename.strip()

        if filename.upper().endswith('.MCV') or filename.upper().endswith('.MC') :
            binary = True
            self.fp = open(filename,"rb")
        else:
            binary = False
            self.fp = open(filename,"r")
        fp= self.fp

        self.name=os.path.splitext(filename)[0]
        # read curve version (INTEGER)
        if binary:
            str = self.getString(4) # dummy 4 '\4\0\0\0'
            self.version_mc_curve = self.getInteger()
            str = self.getString(8) # dummy 8 '\4\0\0\0' + '(' + '\0\0\0'
        else:
            self.version_mc_curve = int(self.fp.readline().strip())
    
        # read dummy text and title 2x (CHARACTER*40)
        if binary:
            str = self.getString(40) # info text '*** File with magnetic curve ***'
            str = self.getString(8) # dummy 8 '(('
            self.mc1_title = self.getString(40) # read title
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
                    energy.append( f )
            else:
                iLen = min(int(s) for s in [self.MC1_MIMAX, (len(values)-idxOffset)])
                for I in range(iLen):
                    energy.append( values[idxOffset+I] )
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

if __name__ == "__main__":
    if len(sys.argv)==2 :
        filename=sys.argv[1]
    else:
        filename=sys.stdin.readline().strip()

    mcv = McvReader()
    mcv.readMcv(filename)
    json.dump( mcv.get_results(), sys.stdout )
