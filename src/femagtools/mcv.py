# -*- coding: utf-8 -*-
"""Reading, Creating and managing MCV/MC files

"""
import json
import functools
import sys
import copy
import logging
import os.path
import pathlib
import struct
import math
import numpy as np
from scipy.interpolate import make_interp_spline
from six import string_types
import femagtools.losscoeffs as lc

# curve types
types = {1: 'Soft iron B(H)',
         2: 'Permanent magnet B(H)',
         3: 'Soft iron B(H,alfa)',
         4: 'Permanent magnet B(H,Br)',
         5: 'Permanent magnet B(H,alfa)',
         6: 'Soft iron Punching',
         7: 'Soft iron Tension',
         8: 'Soft iron Temperature'}

MAGCRV = 1
DEMCRV = 2
ORIENT_CRV = 3
DEMCRV_BR = 4
ORIENT_PM_CRV = 5
PUNCH_CRV = 6
TENSION_CRV = 7
TEMP_CRV = 8
MAG_AC_CRV = -1

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
    ce='mc1_ce_factor',
    ch_freq='mc1_ch_freq_factor',
    cw_freq='mc1_cw_freq_factor',
    fillfac='mc1_fillfac',
    fillfac_old='mc1_fillfac_old',
    bref='mc1_bref',
    bsat='mc1_bsat',
    Bo='mc1_base_induction',
    b_coeff='mc1_induction_factor',
    b_beta_coeff='mc1_induction_beta_factor',
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


def norm_pfe(B, pfe):
    """normalize B and pfe
    B. list of list of flux induction values
    pfe: list of list of loss values
    returns B vector and pfe matrix
    """
    bmin = np.ceil(10*max([b[0] for b in B]))/10.0
    bmax = round(10*max([b[-1] for b in B]))/10.0
    Bv = np.arange(bmin, bmax+0.01, 0.1)
    m = []
    for i, b in enumerate(B):
        n = len([x for x in Bv if x < b[-1]])
        if n < len(b) and n < len(Bv):
            if Bv[n] < b[-1]+0.01:
                b = list(b)
                b[-1] = Bv[n]
                n += 1
        k = 3 if len(pfe[i]) > 3 else 2
        pfunc = make_interp_spline(b, pfe[i], k=k)
        m.append([float(pfunc(x))
                  for x in Bv[:n]] + [None]*(len(Bv)-n))
    return Bv.tolist(), m

def extrapol(crv, jsat=0, bmax=3):
    """extends BH curve into saturation range"""
    curve = copy.deepcopy(crv)
    dh = curve['hi'][-1]-curve['hi'][-2]
    db = curve['bi'][-1]-curve['bi'][-2]
    mue_d = db/dh
    mue = curve['bi'][-1]/curve['hi'][-1]
    db = 3e-2*curve['bi'][-1]
    b2 = 1.5

    curve['muer'] = [b/MUE0/h if h>0 else float('nan')
                     for b, h in zip(curve['bi'],
                                     curve['hi'])]

    # extend bh-curve into saturation
    while mue_d > 1.01*MUE0 and mue > 1.5*MUE0:
        mue_d = MUE0 + (mue_d-MUE0)/np.sqrt(b2)
        curve['bi'].append(curve['bi'][-1]+db)
        curve['hi'].append(curve['hi'][-1]+db/mue_d)
        curve['muer'].append(curve['bi'][-1]/MUE0/curve['hi'][-1])
        b2 += 2
        mue = curve['bi'][-1]/curve['hi'][-1]

    # Fröhlich-Kennelly coefficients
    if jsat:
        hx = curve['hi'][-1]
        bx = curve['bi'][-1]
        fkb = 1./jsat
        fka = (hx/(bx - MUE0*hx) - hx/jsat)
    else:
        Js1 = curve['bi'][-1] - MUE0*curve['hi'][-1]
        Js0 = curve['bi'][-2] - MUE0*curve['hi'][-2]
        fkb = ((curve['hi'][-1]/Js1 - curve['hi'][-2]/Js0)
               /(curve['hi'][-1] - curve['hi'][-2]))
        fka = curve['hi'][-1]/Js1 - fkb*curve['hi'][-1]
    logger.debug("B values %d Fröhlich-Kennelly coeffs: a %f b %f",
                 len(curve['bi']), fka, fkb)
    bstep = 0.15
    bx = np.arange(curve['bi'][-1] + bstep, bmax, bstep)
    # Fröhlich-Kennelly approximation
    b = fkb * bx - MUE0*fka - 1
    a = fkb * bx
    c = MUE0*fka
    nuer = (b + np.sqrt(b*b + 4*a*c))/2/a
    curve['bi'] += bx.tolist()
    curve['hi'] += (nuer*bx/MUE0).tolist()
    return curve


def recalc_bsin(curve):
    """recalculates B-H curve (static problems) into effective
       B-H curve for eddy current problems (voltage driven)."""
    ncurve = []
    ndel = 80
    x = np.linspace(1, ndel, ndel)
    for c in curve:
        nc = dict(
            bi=c['bi'],
            hi=c['hi'][:2])
        bi = list(c['bi'])
        hi = list(c['hi'])
        if bi[0] > 0:
            bi.insert(0, 0)
            hi.insert(0, 0)
        bh = make_interp_spline(bi, hi)
        for bx in c['bi'][2:]:
            bt = bx*np.sin(2*np.pi/4/ndel*x)
            nue = np.sum(bh(bt)/bt)/ndel
            nc['hi'].append(bx*nue)
        # nc['hi'].append(c['hi'][-1])
        ncurve.append(nc)
    return ncurve


def recalc_hsin(curve):
    """recalculates B-H curve (static problems) into effective
       B-H curve for eddy current problems (current driven)."""
    ncurve = []
    ndel = 80
    Tp = 1/200
    dt = Tp/ndel
    x = np.linspace(1, ndel, ndel)
    for c in curve:
        nc = dict(
            hi=c['hi'],
            bi=c['bi'][:2])
        hi = list(c['hi'])
        bi = list(c['bi'])
        if hi[0] > 0:
            hi.insert(0, 0)
            bi.insert(0, 0)
        hb = make_interp_spline(hi, bi)
        for hx in c['hi'][2:]:
            ht = hx*np.sin(2*np.pi/4/ndel*x)
            bt = hb(ht)*np.sin(2*np.pi/4/ndel*x)
            nc['bi'].append(
                2*np.sum((bt[:-1] + bt[1:])*dt/2)/Tp)
        # nc['bi'].append(c['bi'][-1])
        ncurve.append(nc)
    return ncurve


def approx(db2, curve, ctype):
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

    if ctype in (DEMCRV, MAG_AC_CRV):
        dhdbn = 0
        k = len(bi2)-1
        if k < len(curve['bi']):
            if curve['bi'][k] - curve['bi'][k-1] > 0:
                dhdbn = ((curve['hi'][k] - curve['h'][k-1])
                         /(curve['bi'][k] - curve['bi'][k-1]))
            a.append(MUE0*dhdbn)
            b.append(MUE0*curve['hi'][k] - dhdbn*curve['bi'][k])
    else:
        a.append(1.0)
        b.append(MUE0*curve['hi'][-1]-curve['bi'][-1])
    return dict(nuer=nuer, a=a, b=b, bi2=bi2)


def fe_sat_mag(curve):
    """returns maximum polarization of all BH curves"""
    fesat = 0
    for c in curve:
        js2 = c['bi'][-1] - MUE0*c['hi'][-1]
        js1 = c['bi'][-2] - MUE0*c['hi'][-2]
        s = (js1 + js2)/2
        if s > fesat:
            fesat = s
    return fesat


def findNotNone(l):
    """return lower and upper indexes of not none values in list"""
    for i in range(len(l)):
        if l[i]:
            break
    for j in range(len(l)-1, -1, -1):
        if l[j]:
            break
    return (i, j)


class Mcv(object):
    def __init__(self, data={}):
        # default values from file: mcv.par
        self.ACT_VERSION_MC_CURVE = 0
        self.ORIENTED_VERSION_MC_CURVE = 1
        self.PARAMETER_PM_CURVE = 2

        self.MC1_BASE_FREQUENCY = 50.0
        self.MC1_BASE_INDUCTION = 1.5
        self.MC1_CH_FACTOR = 0.0
        self.MC1_CW_FACTOR = 0.0
        self.MC1_CE_FACTOR = 0.0
        self.MC1_CH_FREQ_FACTOR = 0.0
        self.MC1_CW_FREQ_FACTOR = 0.0
        self.MC1_INDUCTION_FACTOR = 0.0
        self.MC1_INDUCTION_BETA_FACTOR = 0.0
        self.jordan = {}
        # {'ch': 0, 'cw': 0, 'ch_freq':0, 'cw_freq':0}
        self.steinmetz = {}
        # {'ch': 0, 'cw': 0, 'ch_freq':0, 'cw_freq':0}
        self.bertotti = {}
        # {'ch': 0, 'cw': 0, 'ce':0, 'ch_freq':0, 'cw_freq':0}
        self.steinmetz_modified = {}
        # {'ch': 0, 'cw': 0, 'ch_freq': 1, 'b_beta_coeff': 0, 'cw_freq': 2, 'b_coeff': 2}
        self.MC1_FE_SPEZ_WEIGTH = 7.65
        self.MC1_FE_SAT_MAGNETIZATION = 2.15

        self.mc1_base_frequency = self.MC1_BASE_FREQUENCY
        self.mc1_base_induction = self.MC1_BASE_INDUCTION
        self.mc1_ch_factor = self.MC1_CH_FACTOR
        self.mc1_cw_factor = self.MC1_CW_FACTOR
        self.mc1_ce_factor = self.MC1_CE_FACTOR
        self.mc1_ch_freq_factor = self.MC1_CH_FREQ_FACTOR
        self.mc1_cw_freq_factor = self.MC1_CW_FREQ_FACTOR
        self.mc1_induction_factor = self.MC1_INDUCTION_FACTOR
        self.mc1_induction_beta_factor = self.MC1_INDUCTION_BETA_FACTOR
        self.mc1_fe_spez_weigth = self.MC1_FE_SPEZ_WEIGTH

        self.mc1_title = ''
        self.version_mc_curve = self.ACT_VERSION_MC_CURVE
        self.mc1_type = MAGCRV    # Soft Iron B(H)
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

        if data:
            self.setData(data)

            self.mc1_curves = len(self.curve)
            if (self.mc1_type in (ORIENT_CRV, ORIENT_PM_CRV)
                or self.mc1_curves > 1):
                self.version_mc_curve = self.ORIENTED_VERSION_MC_CURVE
            elif self.mc1_type == DEMCRV_BR:
                self.version_mc_curve = self.PARAMETER_PM_CURVE


    def __setattr__(self, key, val):
        try:
            self.__dict__[key] = val
        except Exception:  # file format unknown
            logger.debug("setAttr Exception, name: %s, value: %s", key, val)

    def setData(self, data):
        wtrans = {transl[k]: k
                  for k in transl if not isinstance(transl[k], list)}
        for k in wtrans:
            if wtrans[k] in data.keys():
                self.__setattr__(k, data[wtrans[k]])
        for k in ('bertotti', 'jordan', 'steinmetz', 'steinmetz_modified'):
            if k in data:
                self.__setattr__(k, data[k])
        self.curve = data['curve']
        try:
            self.mc1_angle = [c['angle'] for c in data['curve']]
        except Exception:
            self.mc1_angle = [0]*len(data['curve'])
        try:
            self.losses = data['losses']
        except Exception:
            pass
        # assume jordan iron loss parameters
        for k in self.jordan:
            self.jordan[k] = getattr(self, transl[k])
        # set loss coeffs when existing
        try:
            for k in ('bertotti', 'jordan', 'steinmetz', 'steinmetz_modified'):
                if k in data:
                    self.__setattr__(k, data[k])
        except:
            pass
        return

    def rtrimValueList(self, vlist):
        """cut list at first 0"""
        le = len(vlist)
        for i in range(le-1, -1, -1):
            if vlist[i] != 0.:
                break
        return list(vlist[:i+1])

    def __getitem__(self, key):
        if key == 'ctype':  # for compatibility purposes
            return self.mc1_type
        return getattr(self, key)


class Writer(Mcv):
    def __init__(self, data=None):
        Mcv.__init__(self, data)

    def getBlockLength(self, d):
        if isinstance(d, string_types) or isinstance(d, bytes):
            try:
                s = bytes(d).decode('utf-8').encode('latin1')
            except Exception:
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
            except Exception:
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

    def _prepare(self, fillfac, recsin):
        """prepare output format (internal use only)"""
        curve = copy.deepcopy(self.curve)
        if fillfac:
            alpha = fillfac/self.mc1_fillfac
            for c in curve:
                c['bi'] = [alpha*b + MUE0*(1. - alpha)*h
                           for b, h in zip(c['bi'], c['hi'])]
        if recsin == 'flux':
            curve = recalc_bsin(curve)
        elif recsin == 'cur':
            curve = recalc_hsin(curve)
        if fillfac or recsin:
            if hasattr(self, 'mc1_fe_sat_magnetization'):
                self.mc1_fe_sat_magnetization = fe_sat_mag(curve)
        logger.info("%s Type: %d (%s) Num Curves %d",
                    self.name, self.version_mc_curve,
                    types[self.mc1_type],
                    len(self.curve))
        self.mc1_ni = [min(len(c['hi']),
                           len(c['bi']))
                       for c in self.curve if 'hi' in c]
        self.mc1_db2 = [(c['bi'][-1]**2 - c['bi'][0]**2)/n
                        for c, n in zip(curve, self.mc1_mi)]
        for db2, c in zip(self.mc1_db2, curve):
            c.update(approx(db2, c, self.mc1_type))
        if not hasattr(self, 'mc1_fe_sat_magnetization'):
            self.mc1_fe_sat_magnetization = fe_sat_mag(curve)
        self.mc1_mi = [len(c['a'])
                       for c in curve]
        return curve

    def writeBinaryFile(self, fillfac=None, recsin='', feloss=''):
        """write binary file after conversion if requested.
        arguments:
        fillfac: (float) fill actor
        recsin: (str) either 'flux' or 'cur'
        feloss: (str) iron loss method (bertotti, jordan)
        """
        curve = self._prepare(fillfac, recsin)
        write_losses = True
        try:
            if feloss.lower() == 'bertotti':
                for k in self.bertotti:
                    setattr(self, transl[k], self.bertotti[k])
                write_losses = False
            elif feloss.lower() == 'jordan':
                for k in self.jordan:
                    setattr(self, transl[k], self.jordan[k])
                write_losses = False
            elif feloss.lower() == 'steinmetz':
                for k in self.steinmetz:
                    setattr(self, transl[k], self.steinmetz[k])
                write_losses = False
            elif feloss.lower() == 'modified steinmetz':
                for k in self.steinmetz_modified:
                    setattr(self, transl[k], self.steinmetz_modified[k])
                write_losses = False
        except AttributeError as e:
            logger.warning("%s", e)
            pass
        mc1_type = self.mc1_type
        mc1_recalc = self.mc1_recalc
        mc1_fillfac = self.mc1_fillfac
        if fillfac and fillfac < 1:
            mc1_recalc = 1
            mc1_fillfac = fillfac
        if recsin in ('flux', 'cur'):
            mc1_type = MAG_AC_CRV
            mc1_recalc = 1
        # write line, version_mc_curve
        self.writeBlock(self.version_mc_curve)

        # write line, text '    *** File with magnetic curve ***    '
        self.writeBlock('    *** File with magnetic curve ***    ')

        # write line, mc1_title
        self.writeBlock(self.mc1_title.ljust(40)[:40])
        # write line, mc1_ni(1),mc1_mi(1),mc1_type,mc1_recalc,mc1_db2(1)
        self.writeBlock([int(self.mc1_ni[0]),
                         int(self.mc1_mi[0]),
                         int(mc1_type),
                         int(mc1_recalc),
                         self.mc1_db2[0]])

        # write line, mc1_remz, mc1_bsat, mc1_bref, mc1_fillfac
        if self.version_mc_curve == self.ACT_VERSION_MC_CURVE:
            self.writeBlock([float(self.mc1_remz), float(self.mc1_bsat),
                             float(self.mc1_bref), float(mc1_fillfac)])
        if mc1_type == DEMCRV_BR:
            self.mc1_remz = self.mc1_angle[self.mc1_curves-1]
        if self.version_mc_curve in (self.ORIENTED_VERSION_MC_CURVE,
                                     self.PARAMETER_PM_CURVE):
            logging.debug("write mc1_curves %d", self.mc1_curves)
            self.writeBlock([float(self.mc1_remz), float(self.mc1_bsat),
                             float(self.mc1_bref), float(mc1_fillfac),
                             self.mc1_curves])

        if mc1_type == DEMCRV_BR:
            self.mc1_angle[self.mc1_curves-1] = self.mc1_remz

        # data
        for K in range(0, self.mc1_curves):
            logger.debug(" K %d  Bi, Hi %d", K,  len(curve[K].get('bi', [])))
            # hi, bi
            lb = curve[K].get('bi', [])
            lh = curve[K].get('hi', [])
            self.writeBlock(zip(*[
                [float(lb[j]) if j < len(lb) else 0.
                 for j in range(self.MC1_NIMAX)],
                [float(lh[j]) if j < len(lh) else 0.
                 for j in range(self.MC1_NIMAX)]]))

            # bi2, nuer
            lb = curve[K]['bi2']
            ln = curve[K]['nuer']
            self.writeBlock(zip(*[
                [float(lb[j]) if j < len(lb) else 0.
                 for j in range(self.MC1_NIMAX)],
                [float(ln[j]) if j < len(ln) else 0.
                 for j in range(self.MC1_NIMAX)]]))

            # a, b, c, d
            la = curve[K].get('a', [0.]*self.MC1_NIMAX)
            lb = curve[K].get('b', [0.]*self.MC1_NIMAX)
            self.writeBlock(zip(*[
                [float(la[j]) if j < len(la) else 0.
                 for j in range(self.MC1_NIMAX)],
                [float(lb[j]) if j < len(lb) else 0.
                 for j in range(self.MC1_NIMAX)],
                [0.]*50,
                [0.]*50
            ]))

            if self.version_mc_curve in (self.ORIENTED_VERSION_MC_CURVE,
                                         self.PARAMETER_PM_CURVE):
                self.writeBlock([self.mc1_angle[K], self.mc1_db2[K]])

        try:
            if (not (self.mc1_ch_factor or self.mc1_cw_factor)
                and self.losses and write_losses):
                # fit loss parameters
                pfe = self.losses['pfe']
                f = self.losses['f']
                B = self.losses['B']
                losses = [list(p) + [None]*(len(B)-len(p)) for p in pfe]
                fo = self.mc1_base_frequency
                Bo = self.mc1_base_induction
                fit_jordan = False
                if fit_jordan:
                    ch, alfa, cw, beta, gamma = lc.fitjordan(f, B, losses, Bo, fo)
                    self.mc1_ch_factor = ch
                    self.mc1_ch_freq_factor = alfa
                else:
                    cw, beta, gamma = lc.fitsteinmetz(f, B, losses, Bo, fo)
                    self.mc1_ch_factor = 0
                    self.mc1_ch_freq_factor = 0
                self.mc1_cw_factor = cw
                self.mc1_cw_freq_factor = beta
                self.mc1_induction_factor = gamma

        except AttributeError:
            pass
        self.writeBlock([float(self.mc1_base_frequency),
                         float(self.mc1_base_induction),
                         float(self.mc1_ch_factor),
                         float(self.mc1_cw_factor),
                         float(self.mc1_ch_freq_factor),
                         float(self.mc1_cw_freq_factor),
                         float(self.mc1_induction_factor),
                         float(self.mc1_fe_spez_weigth),
                         float(self.mc1_fe_sat_magnetization)])

        logger.info("MCV coeffs %s",
                    {"fo": float(self.mc1_base_frequency),
                     "Bo": float(self.mc1_base_induction),
                     "ch": float(self.mc1_ch_factor),
                     "cw": float(self.mc1_cw_factor),
                     "ch_freq": float(self.mc1_ch_freq_factor),
                     "cw_freq": float(self.mc1_cw_freq_factor),
                     "b_coeff": float(self.mc1_induction_factor),
                     "fr_spez_weight": float(self.mc1_fe_spez_weigth),
                     "fe_sat_magnetization": float(self.mc1_fe_sat_magnetization)})

        if not hasattr(self, 'losses') or not self.losses:
            # new variables: ce factor for bertotti losses
            # b_beta_coeff for modified steinmetz
            try:
                self.writeBlock([float(self.mc1_ce_factor),
                                 float(self.mc1_induction_beta_factor)])
                #logger.info(f"ce = {float(self.mc1_ce_factor)}")
                #logger.info(f"b_beta_coeff = {float(self.mc1_induction_beta_factor)}")
            except:
                pass
            return

        try:
            freq = [x for x in self.losses['f'] if x > 0]
            nfreq = len(freq)
            nind = len(self.losses['B'])
            if nind < 1 or nfreq < 1:
                return
            fo = self.losses.get('fo',
                                 self.mc1_base_frequency)
            Bo = self.losses.get('Bo',
                                 self.mc1_base_induction)
            if np.isscalar(self.losses['B'][0]):
                B = self.losses['B']
                pfe = self.losses['pfe']
                if 'cw' not in self.losses:
                    cw, alfa, beta = lc.fitsteinmetz(
                        self.losses['f'], self.losses['B'], self.losses['pfe'], Bo, fo)
                    self.losses['cw'] = cw
                    self.losses['cw_freq'] = alfa
                    self.losses['b_coeff'] = beta
                    self.losses['Bo'] = Bo
                    self.losses['fo'] = fo
            else:
                cw, alfa, beta = lc.fitsteinmetz(
                    self.losses['f'], self.losses['B'], self.losses['pfe'], Bo, fo)
                B, pfe = norm_pfe(self.losses['B'], self.losses['pfe'])
                nind = len(B)
                self.losses['cw'] = cw
                self.losses['cw_freq'] = alfa
                self.losses['b_coeff'] = beta
                self.losses['Bo'] = Bo
                self.losses['fo'] = fo

            self.writeBlock([nfreq, nind])
            self.writeBlock([float(b) for b in B] + [0.0]*(M_LOSS_INDUCT - nind))

            nrec = 0
            for f, p in zip(self.losses['f'], pfe):
                if f > 0:
                    y = np.array(p)
                    losses = [float(x) for x in y[y != np.array(None)]]
                    nloss = len(losses)
                    if nloss == nind:
                        pl = list(p)
                    else:
                        cw, alfa, beta = lc.fitsteinmetz(
                            f, B[:nloss], losses, Bo, fo)
                        pl = losses + [lc.pfe_steinmetz(
                            f, b, cw, alfa, beta,
                            self.losses['fo'],
                            self.losses['Bo'])
                                       for b in B[nloss:]]
                    logger.debug("%s", pl)
                    self.writeBlock(pl +
                                    [0.0]*(M_LOSS_INDUCT - len(pl)))
                    self.writeBlock(float(f))
                    nrec += 1

            if write_losses:
                logger.info("Append empty blocks %d x %d",
                            M_LOSS_FREQ - nrec, M_LOSS_INDUCT)
                for m in range(M_LOSS_FREQ - nrec):
                    self.writeBlock([0.0]*M_LOSS_INDUCT)
                    self.writeBlock(0.0)

            self.writeBlock([self.losses['cw'], self.losses['cw_freq'],
                             self.losses['b_coeff'], self.losses['fo'],
                             self.losses['Bo']])
            logger.info("loss coeff %s",
                        {"cw": float(self.losses['cw']),
                         "cw_freq": float(self.losses['cw_freq']),
                         "b_coeff": float(self.losses['b_coeff']),
                         "fo": float(self.losses['fo']),
                         "Bo": float(self.losses['Bo'])})

            self.writeBlock([1])
            logger.info('Losses n freq %d n ind %d', nfreq, nind)
        except Exception as e:
            logger.error("Exception %s", e, exc_info=True)

    def writeMcv(self, filename, fillfac=None, recsin='', feloss='Jordan'):
        # windows needs this strip to remove '\r'
        filename = pathlib.Path(filename)
        self.name = filename.stem

        if filename.suffix.upper() in ('.MCV', '.MC'):
            binary = True
            self.fp = filename.open(mode="wb")
        else:
            binary = False
            self.fp = filename.open(mode="w")
        logger.info("Write File %s, binary format (feloss '%s')",
                    filename, feloss)

        self.writeBinaryFile(fillfac, recsin, feloss)
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
        filename = pathlib.Path(filename)

        if filename.suffix in ('.MCV', '.MC'):
            binary = True
            self.fp = filename.open(mode="rb")
        else:
            binary = False
            self.fp = filename.open(mode="r")

        self.name = filename.stem
        # read curve version (INTEGER)
        if binary:
            self.version_mc_curve = self.readBlock(int)
        else:
            self.version_mc_curve = int(self.fp.readline().strip())

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
            self.mc1_ni[0] = int(line[0])   # mc1_ni (INTEGER)
            self.mc1_mi[0] = int(line[1])   # mc1_mi (INTEGER)
            self.mc1_type = int(line[2])   # mc1_type (INTEGER)
            self.mc1_recalc = int(line[3])   # mc1_recalc (INTEGER)
            self.mc1_db2[0] = float(line[4])  # mc1_db2 (REAL)

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
        logger.info("Read file %s %s Type: %d (%s) Num Curves %d",
                    filename, self.name, self.version_mc_curve,
                    types[self.mc1_type], self.mc1_curves)
        if self.mc1_type == DEMCRV_BR:
            self.mc1_angle[self.mc1_curves-1] = self.mc1_remz

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
                hi=self.rtrimValueList(mc_hi[:self.MC1_NIMAX]),
                bi=self.rtrimValueList(mc_bi[:self.MC1_NIMAX]),
                bi2=self.rtrimValueList(mc_bi2[:self.MC1_MIMAX]),
                nuer=self.rtrimValueList(mc_nuer[:self.MC1_MIMAX]),
                a=self.rtrimValueList(mc_a[:self.MC1_MIMAX]),
                b=self.rtrimValueList(mc_b[:self.MC1_MIMAX])))

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

        if vals[0]:
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

        # TODO: handle self.mc1_ce_factor, self.mc1_induction_beta_factor

        self.losses = {}
        try:
            (nfreq, njind) = self.readBlock([int, int])
            if (nfreq and njind):
                self.losses['B'] = self.readBlock(
                    [float]*M_LOSS_INDUCT)[:njind]
                self.losses['f'] = []
                self.losses['pfe'] = []
                for i in range(M_LOSS_FREQ):
                    res = self.readBlock([float]*M_LOSS_INDUCT)
                    f = self.readBlock(float)
                    if i<nfreq and f != None:
                        self.losses['pfe'].append(res[:njind])
                        self.losses['f'].append(f)
                    else:
                        self.losses['f'].append(0)
                        self.losses['pfe'].append([0]*njind)
                (cw, alfa, beta, basefreq, baseind) = self.readBlock([float]*5)
                self.losses['fo'] = basefreq
                self.losses['Bo'] = baseind
                self.losses['cw'] = cw
                self.losses['cw_freq'] = alfa
                self.losses['b_coeff'] = beta
                self.losses['ch'] = self.ch
                self.losses['ch_freq'] = self.ch_freq
        except Exception as e:
            logger.debug("Exception %s", e)
            if self.losses and 'B' in self.losses:
                if not self.losses['f'] or not self.losses['pfe']:
                    self.losses = {}

        self.ce = 0
        self.b_beta_coeff = 0
        try:
            if not self.losses['f'][0]:
                self.fp.seek(-16, 1)
            res = self.readBlock([float]*2)
            self.ce, self.b_beta_coeff = res[0], res[1]
        except:
            pass

    def get_results(self):
        self.desc = self.mc1_title
        self.fillfac = self.mc1_fillfac
        #self.cversion = self.version_mc_curve
        result = {
            k: getattr(self, k) for k in (
                'name',
                'desc',
                'cversion',
                'ctype',
                'recalc',
                'remz',
                'bsat',
                'bref',
                'fillfac',
                'fo',
                'Bo',
                'ch',
                'ch_freq',
                'cw',
                'ce',
                'b_beta_coeff',
                'cw_freq',
                'b_coeff',
                'rho',
                'fe_sat_mag') if hasattr(self, k)}
        result['curve'] = [{k: c[k] for k in ('hi', 'bi')}
                           for c in self.curve]
        try:
            if self.losses:
                result['losses'] = self.losses
        except Exception:
            pass
        if (self.ORIENTED_VERSION_MC_CURVE or self.PARAMETER_PM_CURVE):
            for i in range(len(self.curve)):
                result['curve'][i]['angle'] = self.mc1_angle[i]

        return result

    def keys(self):
        return [k for k in dir(self) if not k.startswith('_')]


class MagnetizingCurve(object):
    """Magnetizing curve
        Args:
          mcvpar: a list of mcv objects or a single mcv or a directory
    """
    def __init__(self, mcvpar):
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
            try:
                self.mcv[str(mcvpar['id'])] = mcvpar
                return
            except KeyError:
                pass
            try:
                self.mcv[mcvpar['name']] = mcvpar
                return
            except KeyError:
                pass

            self.mcv['0'] = mcvpar
            return

        elif isinstance(mcvpar, Reader):
            self.mcv[mcvpar['name']] = mcvpar.get_results()
            return

        elif isinstance(mcvpar, string_types):
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
                logger.info("search file %s in %s", filename,
                            self.mcdirectory)
                if os.access(os.path.join(self.mcdirectory,
                                          filename), os.R_OK):
                    return id
            except AttributeError as ex:
                #logger.warn("Exception %s", ex)
                pass  # no mcdirectory

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
        try:
            return self.mcv[name]
        except Exception:
            pass
        return None

    def fix_name(self, name, fillfac=1.0):
        """return os compatible mcv name including fillfac"""
        if not self.find_by_name(name):
            if fillfac and fillfac < 1.0:
                return "{0}-{1:d}".format(name, int(100*fillfac))
            return name
        repls = {' ': '_', '(': '_', ')': '_', ',': '_'}
        if fillfac and fillfac < 1.0:
            return "{0}-{1:d}".format(
                functools.reduce(lambda a, kv: a.replace(*kv),
                                 repls.items(), name),
                int(100*fillfac))
        return functools.reduce(lambda a, kv: a.replace(*kv),
                                repls.items(), name)

    def writefile(self, name, directory='.',
                  fillfac=None, recsin='', feloss=''):
        """find magnetic curve by name or id and write binary file
        Arguments:
          name: key of mcv dict (name or id)
          directory: destination directory (must be writable)
          fillfac: (float) new fill factor (curves will be recalulated if not None or 0)
          recsin: (str) either 'flux' or 'cur' recalculates for eddy current calculation (dynamic simulation)
          feloss: (str) iron loss calc method ('jordan', 'bertotti', 'steinmetz')

        returns filename if found else None
        """
        ext = '.MC' if sys.platform == 'win32' else '.MCV'
        mcv = self.find_by_name(name)
        if not mcv:
            bname = name
            filename = ''.join((name, ext))
            # check fillfac and readmcv
            if not fillfac or fillfac == 1.0:
                try:
                    import shutil
                    logger.info("Copy file %s", filename)
                    shutil.copy(os.path.join(self.mcdirectory,
                                             filename), directory)
                    return filename
                except shutil.SameFileError:
                    return filename
                except Exception:
                    logger.error("MCV %s not found in directory %s",
                                 str(filename), directory)
                    return None
            try:
                mcv = Reader()
                mcv.readMcv(os.path.join(self.mcdirectory,
                                         filename))
            except AttributeError:
                logger.error("MCV %s not found in dict list", name)
                return ''
        bname = self.fix_name(mcv['name'], fillfac)
        filename = ''.join((bname, ext))
        writer = Writer(mcv)
        writer.writeMcv(os.path.join(directory, filename),
                        fillfac=fillfac, recsin=recsin, feloss=feloss)
        return filename

    def fitLossCoeffs(self):
        for m in self.mcv:
            if 'losses' not in self.mcv[m]:
                continue
            losses = self.mcv[m]['losses']
            cw, alfa, beta = lc.fitsteinmetz(
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


class MCVconvert:
    def __init__(self, bhdata):
        self.steinmetz = dict(cw=0, cw_freq=0, b_coeff=0)
        self.jordan = dict(ch=0, ch_freq=0, cw=0, cw_freq=0, b_coeff=0)
        self.bertotti = dict(ch=0, cw=0, ce=0)
        self.steinmetz_modified = dict(ch=0, cw=0, ch_freq=1, b_beta_coeff=0, cw_freq=2, b_coeff=2)
        self.losscalc = None

        B = []
        f = [50.0, 100.0, 200.0, 400.0, 1000.0, 2000.0]
        pfe = []
        jordan = False
        bertotti = False
        modified_steinmetz = False
        flag = False

        if "losses" in bhdata:
            # test if any nan in data
            for i in ("B", "f"):
                if np.any(np.isnan(bhdata["losses"][i])) or \
                np.any(np.isnan(bhdata["losses"][i])):
                    flag = True
            for i in bhdata["losses"]['pfe']:
                if np.any(np.isnan(i)):
                    flag = True

        if 'losses' not in bhdata or flag:
            # check steinmetz or jordan
            bhdata.update({"losses": {"pfe": [], "f":[], "B": []}})
            B = bhdata["curve"][-1]["bi"]
            if "ch" in bhdata:

                if "ce" in bhdata:
                    if bhdata['ce'] > 1e-15:
                        bertotti = True

                    if bhdata["ch"] > 1e-15 and bhdata['ce'] < 1e-15:
                        jordan = True

                else:
                    jordan = True

                #if (conditions for modified steinmetz):
                #   modified_steinmetz = True

            if jordan:
                self.losscalc = 'jordan'
                logger.info("calculating based on jordan...")
                for i in f:
                    pfe.append(lc.pfe_jordan(i, np.array(B), bhdata['ch'], bhdata['ch_freq'], bhdata['cw'],
                                          bhdata['cw_freq'], bhdata['b_coeff'], 50, 1.5))

            elif bertotti:
                self.losscalc = 'bertotti'
                logger.info("calculating based on bertotti")
                for i in f:
                    pfe.append(lc.wbert(i, np.array(B), bhdata['ch'], bhdata['cw'], bhdata['ce']))

            elif modified_steinmetz:
                self.losscalc = 'modified steinmetz'
                logger.info("calculating based on modified steinmetz...")
                for i in f:
                    pfe.append(lc.pfe_modified_steinmetz(i, np.array(B), bhdata['ch'], bhdata['cw'], bhdata['b_coeff'], bhdata['b_beta_coeff'], 1, 1))

            else: # steinmetz
                self.losscalc = 'steinmetz'
                logger.info("calculating based on steinmetz...")
                for i in f:
                    pfe.append(lc.pfe_steinmetz(i, np.array(B), bhdata['cw'],
                                            bhdata['cw_freq'], bhdata['b_coeff'], 50, 1.5))
            bhdata['losses']['pfe'] = pfe
            bhdata['losses']['B'] = B
            bhdata['losses']['f'] = f

        idx = 0
        for i, j in enumerate(bhdata['losses']['f']):
            if j == 0:
                idx = i
                break
        idx = idx - 1
        z = lc.fitsteinmetz(bhdata['losses']['f'][0:idx],
                                     bhdata['losses']['B'],
                                     bhdata['losses']['pfe'][0:idx],
                                     1.5,
                                     50)

        for i, j in enumerate(self.steinmetz):
            self.steinmetz[j] = z[i]
        self.steinmetz.update({"Bo": 1.5, "fo": 50})

        z = lc.fitjordan(bhdata['losses']['f'][0:idx],
                    bhdata['losses']['B'],
                    bhdata['losses']['pfe'][0:idx],
                    1.5,
                    50)

        for i, j in enumerate(self.jordan):
            self.jordan[j] = z[i]
        self.jordan.update({"Bo": 1.5, "fo": 50})

        z = lc.fit_bertotti(bhdata['losses']['f'][0:idx],
                        bhdata['losses']['B'],
                        bhdata['losses']['pfe'][0:idx])

        for i, j in enumerate(self.bertotti):
            self.bertotti[j] = z[i]
        self.bertotti.update({"Bo": 1, "fo": 1, "alpha": 2.0, "ch_freq": 1.0, "cw_freq": 2.0, "b_coeff": 2.0})

        z = lc.fit_modified_steinmetz(bhdata['losses']['f'][0:idx],
                        bhdata['losses']['B'],
                        bhdata['losses']['pfe'][0:idx], 1, 1)

        for i,j in enumerate(self.steinmetz_modified):
            self.steinmetz_modified[j] = z[i]
        self.steinmetz_modified.update({"Bo": 1, "fo": 1})

        bhdata['losses']['pfe'] = np.transpose(bhdata['losses']['pfe']).tolist() #len(B) rows and len(f) columns


def read(filename):
    """read MC/MCV file and return mc dict"""
    mcv = Reader()
    mcv.readMcv(filename)
    return mcv


if __name__ == "__main__":
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    mcv = read(filename)
    json.dump(mcv.get_results(), sys.stdout)
