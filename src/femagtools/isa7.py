# -*- coding: utf-8 -*-
""" Read FEMAG I7/ISA7 model files

"""
import re
import sys
import struct
from enum import Enum
import logging
from collections import Counter, defaultdict
import numpy as np

logger = logging.getLogger('femagtools.isa7')


def Trot(alpha):
    return np.array([[np.cos(alpha), np.sin(alpha)],
                     [-np.sin(alpha), np.cos(alpha)]])

def transform_coord(geometry, xy):
    '''transform from global coord to local coord'''
    ecpl = Trot(geometry['alpha']).dot((xy-geometry['cxy']).T).T
    return dict(ecpl=(ecpl + (geometry['w']/2,
                              geometry['h']/2)).T,
                ecp=np.asarray(xy).T)


def transform_flux_density(alpha, bxy):
    '''transform the magnet flux density to local coordinate system'''
    def tf(b1, b2, alpha):
        if b1.ndim > 1:
            r = Trot(alpha).dot(((b1.ravel()), (b2.ravel())))
            return [r[0, :].reshape(*b1.shape),
                    r[1, :].reshape(*b1.shape)]
        else:
            return Trot(alpha).dot(((b1), (b2)))

    b = tf(b1=bxy[:, 0, :], b2=bxy[:, 1, :], alpha=alpha)

    # remove DC component
    bxf = np.mean(b[0].T - np.mean(b[0], axis=1).T, axis=1)
    byf = np.mean(b[1].T - np.mean(b[1], axis=1).T, axis=1)
    return {'bxyl': np.asarray(b),
            'bxyf': np.array([bxf, byf])}


def jordanpfe(Bxnu, Bynu, fnu, losscoeffs, axr):
    """example of custom core loss calculation
    Args:
        Bxnu, Bynu: float arrays of flux density
        fnu: float arrays of frequency
        losscoeffs: dict with keys
              base_frequency, base_induction,
              ch, cw, ch_freq_exp, cw_freq_exp, cw_ind_exp
    """
    basfrq = losscoeffs['base_frequency']
    basind = losscoeffs['base_induction']
    ch = losscoeffs['ch']
    cw = losscoeffs['cw']
    hyscoef = losscoeffs['ch_freq_exp']
    edycoef = losscoeffs['cw_freq_exp']
    indcoef = losscoeffs['cw_ind_exp']
    br = (np.sqrt(Bxnu**2+Bynu**2)/basind)**indcoef
    fr = fnu/basfrq
    hch, hcw = fr**hyscoef, fr**edycoef
    phy, pec = ch*hch*br, cw*hcw*br
    return phy, pec, (0,)

def bertotti_pfe(Bxnu, Bynu, fnu, losscoeffs, axr):
    """bertotti core loss formula
    """
    basfrq = losscoeffs['base_frequency']
    basind = losscoeffs['base_induction']

    ch = losscoeffs['ch'] # hysteresis
    cw = losscoeffs['cw'] # eddy current
    ce = losscoeffs['ce'] # excess

    '''
    hyscoef = losscoeffs['ch_freq_exp']
    edycoef = losscoeffs['cw_freq_exp']
    exccoef = losscoeffs['cx_freq_exp']
    hys_indcoef = losscoeffs['ch_ind_exp']
    edy_indcoef = losscoeffs['cw_ind_exp']
    exc_indcoef = losscoeffs['cx_ind_exp']
    '''

    b21 = np.linalg.norm((Bxnu, Bynu), axis=0)
    b = (b21/basind)
    hi = fnu/basfrq

    # classic bertotti formula
    phy = ch*hi*b
    pec = cw*(hi**2)*(b**2)
    pex = ce*(hi**1.5)*(b**1.5)

    return phy, pec, pex

class ElType(Enum):
    LinearTriangle = 1
    LinearRectangle = 2
    SquareTriangle = 3
    SquareRectangle = 4
"""Types of elements"""

class MoveType(Enum):
    Rotate = 0
    Linear = 1
"""Types of rotor movement"""

class Reader(object):
    """
    Open and Read I7/ISA7 file

    Arguments:
        filename: name of I7/ISA7 file to be read
    """

    def __init__(self, filename):
        self.BR_TEMP_COEF = 0
        with open(filename, mode="rb") as self.file:
            self.file = self.file.read()
        self.pos = 0
        (self.NUM_PNT, self.PNT_PTR, self.PNT_HIDX,
         self.NUM_LIN, self.LIN_PTR, self.LIN_HIDX,
         self.NUM_NOD, self.NOD_PTR, self.NOD_HIDX,
         self.NUM_NDEL, self.NDEL_PTR, self.NDEL_HIDX,
         self.NUM_NDCH, self.NDCH_PTR, self.NDCH_HIDX,
         self.NUM_ELE, self.ELE_PTR, self.ELE_HIDX,
         self.NUM_ELND, self.ELND_PTR, self.ELND_HIDX,
         self.NUM_SPEL, self.SPEL_PTR, self.SPEL_HIDX,
         self.NUM_SE_EL, self.SE_EL_PTR, self.SE_EL_HIDX,
         self.NUM_SPEL_NDCH, self.SPEL_NDCH_PTR, self.SPEL_NDCH_HIDX,
         self.NUM_SR, self.SR_PTR, self.SR_HIDX,
         self.NUM_SR_SE, self.SR_SE_PTR, self.SR_SE_HIDX,
         self.NUM_WB, self.WB_PTR, self.WB_HIDX,
         self.NUM_WB_SR, self.WB_SR_PTR, self.WB_SR_HIDX,
         self.NUM_OB, self.OB_PTR, self.OB_HIDX,
         self.NUM_OB_SR, self.OB_SR_PTR, self.OB_SR_HIDX,
         self.NUM_DV, self.DV_PTR, self.DV_HIDX,
         self.NUM_DV_OB, self.DV_OB_PTR, self.DV_OB_HIDX,
         self.NUM_MC, self.MC_PTR, self.MC_HIDX,
         self.NUM_CF, self.CF_PTR, self.CF_HIDX,
         self.NUM_CF_MC, self.CF_MC_PTR, self.CF_MC_HIDX,
         self.NUM_WN, self.WN_PTR, self.WN_HIDX,
         self.NUM_WN_SW, self.WN_SW_PTR, self.WN_SW_HIDX
         ) = self.next_block("i")

        (valid,
         self.POINT_ISA_POINT_REC_PT_CO_X,
         self.POINT_ISA_POINT_REC_PT_CO_Y) = self.next_block("?ff")

        (valid,
         self.LINE_ISA_LINE_REC_LN_PNT_1,
         self.LINE_ISA_LINE_REC_LN_PNT_2) = self.next_block("?hh")

        (valid,
         self.NODE_ISA_NOD_EL_PNTR,
         self.NODE_ISA_ND_CO_RAD,
         self.NODE_ISA_ND_CO_PHI,
         self.NODE_ISA_NODE_REC_ND_BND_CND,
         self.NODE_ISA_NODE_REC_ND_PER_NOD,
         self.NODE_ISA_NODE_REC_ND_SV_PNTR,
         self.NODE_ISA_NODE_REC_ND_CO_1,
         self.NODE_ISA_NODE_REC_ND_CO_2,
         self.NODE_ISA_NODE_REC_ND_VP_RE,
         self.NODE_ISA_NODE_REC_ND_VP_IM) = self.next_block("?iffhiiffff")

        (self.NOD_ELE_ISA_EL_KEY,
         self.NOD_ELE_ISA_NXT_EL_PNTR) = self.next_block("ii")

        (valid,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_1,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_2,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_MID) = self.next_block("?iii")

        (valid,
         self.ELEM_ISA_EL_NOD_PNTR,
         self.ELEM_ISA_ELEM_REC_EL_TYP,
         self.ELEM_ISA_ELEM_REC_EL_SE_KEY,
         self.ELEM_ISA_ELEM_REC_EL_RELUC,
         self.ELEM_ISA_ELEM_REC_EL_RELUC_2,
         self.ELEM_ISA_ELEM_REC_EL_MAG_1,
         self.ELEM_ISA_ELEM_REC_EL_MAG_2) = self.next_block("?ihhffff")

        (self.ELE_NOD_ISA_ND_KEY,
         self.ELE_NOD_ISA_NXT_ND_PNTR) = self.next_block("ii")

        (valid,
         self.SUPEL_ISA_SE_NDCHN_PNTR,
         self.SUPEL_ISA_SE_EL_PNTR,
         self.SUPEL_ISA_SUPEL_REC_SE_COL,
         self.SUPEL_ISA_SUPEL_REC_SE_MCV_TYP,
         self.SUPEL_ISA_SUPEL_REC_SE_COND_TYP,
         self.SUPEL_ISA_SUPEL_REC_SE_VEL_SYS,
         self.SUPEL_ISA_SUPEL_REC_SE_SR_KEY,
         self.SUPEL_ISA_SUPEL_REC_SE_VELO_1,
         self.SUPEL_ISA_SUPEL_REC_SE_VELO_2,
         self.SUPEL_ISA_SUPEL_REC_SE_CONDUC,
         self.SUPEL_ISA_SUPEL_REC_SE_LENGHT,
         self.SUPEL_ISA_SUPEL_REC_SE_CURD_RE,
         self.SUPEL_ISA_SUPEL_REC_SE_CURD_IM
         ) = self.next_block("?iihhhhhffffff")

        (self.SE_NDCHN_ISA_NC_KEY,
         self.SE_NDCHN_ISA_NXT_NC_PNTR) = self.next_block("ii")

        (self.SE_EL_ISA_EL_KEY,
         self.SE_EL_ISA_NXT_EL_PNTR) = self.next_block("ii")

        (valid,
         self.SR_ISA_SR_SE_PNTR,
         self.SR_ISA_SR_REC_SR_TYP,
         self.SR_ISA_SR_REC_SR_COL,
         self.SR_ISA_SR_REC_SR_NAME,
         self.SR_ISA_SR_REC_SR_CUR_DIR,
         self.SR_ISA_SR_REC_SR_WB_KEY,
         self.SR_ISA_SR_REC_SR_NTURNS,
         self.SR_ISA_SR_REC_SR_SV_PNTR,
         self.SR_ISA_SR_REC_SR_ARRAY,
         self.SR_ISA_SR_REC_SR_GCUR_RE,
         self.SR_ISA_SR_REC_SR_GCUR_IM,
         self.SR_ISA_SR_REC_SR_VOLT_RE,
         self.SR_ISA_SR_REC_SR_VOLT_IM) = self.next_block("?hhh4shhhhfffff")

        (self.SR_SE_ISA_SE_KEY,
         self.SR_SE_ISA_NXT_SE_PNTR) = self.next_block("hh")

        (valid,
         self.WB_ISA_WB_SR_PNTR,
         self.WB_ISA_WB_REC_WB_COL,
         self.WB_ISA_WB_REC_WB_NAME,
         self.WB_TURN,
         self.WB_ISA_WB_REC_WB_SR_NUM,
         self.WB_ISA_WB_REC_WB_WND_KEY,
         self.WB_ISA_WB_REC_WB_UNIT_RES,
         self.WB_ISA_WB_REC_WB_GCUR_RE,
         self.WB_ISA_WB_REC_WB_GCUR_IM,
         self.WB_ISA_WB_REC_WB_VOLT_RE,
         self.WB_ISA_WB_REC_WB_VOLT_IM,
         self.WB_ISA_WB_REC_WB_IMPDZ_RE,
         self.WB_ISA_WB_REC_WB_IMPDZ_IM) = self.next_block("?hh4shhhfffffff")

        self.WB_ISA_WB_REC_WB_TURN = []
        for wd in range(self.NUM_WB):
            if self.WB_ISA_WB_REC_WB_UNIT_RES[wd] == 0:
                self.WB_ISA_WB_REC_WB_TURN.append(
                    self.WB_TURN[wd])
            else:
                self.WB_ISA_WB_REC_WB_TURN.append(
                    self.WB_ISA_WB_REC_WB_UNIT_RES[wd])
                self.WB_ISA_WB_REC_WB_UNIT_RES[wd] = 0

        (self.WB_SR_ISA_SR_KEY,
         self.WB_SR_ISA_NXT_SR_PNTR) = self.next_block("hh")

        self.skip_block(21)

        ANZAHL_TG = self.next_block("i")[1]

        self.skip_block(7)
        self.skip_block(ANZAHL_TG + 1)
        self.skip_block(1)

        self.FC_RADIUS = self.next_block("f")[0]
        self.skip_block(2)
        self.M_POLES = self.next_block("i")[0]
        self.skip_block(1)
        length_cu = [0,0]
        fillfactor_cu = [0,0]
        sigma_cu = [0,0]
        length_cu[0], fillfactor_cu[0], sigma_cu[0] = self.next_block("f")[0:3]
        self.skip_block(3)
        self.MAGN_TEMPERATURE, self.BR_TEMP_COEF = self.next_block("f")[0:2]
        FC_NUM_CUR_ID, FC_NUM_BETA_ID, FC_NUM_CUR_REG = self.next_block("i")[0:3]
        if FC_NUM_CUR_ID > 16:
            FC_NUM_CUR_ID = 16

        self.CURRENT_ID = self.next_block("f")[0:FC_NUM_CUR_ID]
        self.IDE_FLUX = self.next_block("f")[0:FC_NUM_CUR_ID]
        self.IDE_BETA = self.next_block("f")[0:FC_NUM_BETA_ID]

        self.skip_block(FC_NUM_CUR_ID * 2)
        self.skip_block(1 + 10 * 5)
        self.HC_TEMP_COEF, self.NHARM_MAX_MULTI_FILE = self.next_block("f")[0:2]
        self.CU_SPEZ_WEIGHT, self.MA_SPEZ_WEIGHT = self.next_block("f")[0:2]
        self.PS_IND_LOW, self.PS_IND_HIGH = self.next_block("f")[0:2]
        self.skip_block(5 + 14)
        NUM_FE_EVAL_MOVE_STEP = self.next_block("i")[0]
        if NUM_FE_EVAL_MOVE_STEP < 0:
            NUM_FE_EVAL_MOVE_STEP = 0

        self.el_fe_induction_1 = []
        self.el_fe_induction_2 = []
        self.eddy_cu_vpot = []
        self.pos_el_fe_induction = []

        if NUM_FE_EVAL_MOVE_STEP > 1:
            self.pos_el_fe_induction = self.next_block("f")
            self.el_fe_induction_1.append([[]])
            self.el_fe_induction_2.append([[]])
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                b = self.next_block("h")
                logger.debug("el_fe_induction move step %d: %d", i, len(b))
                self.el_fe_induction_1[0][0].append(b)
                self.el_fe_induction_2[0][0].append(self.next_block("h"))

        FC_NUM_MOVE_CALC_LOAD_PMS, FC_NUM_FLX = self.next_block("i")[0:2]

        if FC_NUM_MOVE_CALC_LOAD_PMS > 1:
            self.skip_block(4)
            self.skip_block(3 * FC_NUM_FLX)
            self.skip_block()

        FC_NUM_MOVE_NOLOAD_PMS = self.next_block("i")[0]

        if FC_NUM_MOVE_NOLOAD_PMS > 1:
            self.skip_block(4)
            self.skip_block(2 * FC_NUM_FLX)
            self.skip_block()

        if NUM_FE_EVAL_MOVE_STEP > 1:
            self.eddy_cu_vpot.append([[]])
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.eddy_cu_vpot[0][0].append(self.next_block("h"))

        self.skip_block(2)  # start_winkel, end_winkel
        self.skip_block(2 * 5 + 2)
        length_cu[1], fillfactor_cu[1], sigma_cu[1] = self.next_block("f")[0:3]
        self.PS_LENGTH_CU = length_cu
        self.PS_FILFACTOR_CU = fillfactor_cu
        self.PS_SIGMA_CU = sigma_cu

        self.skip_block(12)
        self.skip_block(3 * 30 * 30)
        self.skip_block(3)
        self.skip_block(30 * 30)
        self.skip_block(4)
        # stator 3
        self.skip_block(4)
        (yoke_diam, inside_diam,
         slot_height, slot_h1, slot_h2,
         slot_width, slot_r1, slot_r2) = self.next_block("f")[:8]
        self.skip_block(3)
        # magnet sector
        magn_rad, yoke_rad, magn_height = self.next_block("f")[:3]
        self.da2 = 2*magn_rad*1e-3
        self.dy2 = 2*yoke_rad*1e-3
        self.da1 = inside_diam
        self.dy1 = yoke_diam
        self.skip_block(3)
        # windings generation
        (tot_num_slot, num_phases, num_layers,
         self.NUM_WIRES, self.CURRENT,
         coil_span, num_slots) = self.next_block("f")[:7]
        self.slots = int(tot_num_slot)
        self.num_phases = int(num_phases)
        self.layers = int(num_layers)
        self.slots_gen = int(num_slots)
        self.coil_span = coil_span

        self.skip_block(1)
        (move_action, arm_length, self.SKEW_ANGLE,
         HI, num_move_ar, self.ANGL_I_UP,
         num_par_wdgs, cur_control) = self.next_block("f")[:8]
        self.NUM_PAR_WDGS = int(num_par_wdgs)
        self.move_action = move_action  # rotate 0, linear 1
        self.arm_length = arm_length  # unit is m
        self.skip_block(2)
        self.skip_block(30 * 30)
        self.skip_block(30 * 30)
        self.skip_block(1 * 20)
        self.skip_block(8)
        self.beta_loss = self.next_block(
            "h")[:FC_NUM_BETA_ID]  # BETA_LOSS_EVAL_STEP
        self.curr_loss = self.next_block(
            "h")[:FC_NUM_CUR_ID]  # CURR_LOSS_EVAL_STEP
        FC_NUM_MOVE_LOSSES = self.next_block("i")[0]
        if FC_NUM_MOVE_LOSSES > 1 and NUM_FE_EVAL_MOVE_STEP > 1:
            self.el_fe_induction_1.append([[]])
            self.el_fe_induction_2.append([[]])
            self.eddy_cu_vpot.append([[]])
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                b = self.next_block("h")
                logger.debug("el_fe_induction move losses step %d: %d", i, len(b))
                self.el_fe_induction_1[1][0].append(b)
                self.el_fe_induction_2[1][0].append(self.next_block("h"))
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.eddy_cu_vpot[1][0].append(self.next_block("h"))

        # VIRGIN_PM_SYN
        self.skip_block(3)
        # magnet iron 4
        self.skip_block(3)
        # stator 4
        self.skip_block(1)
        # stator 2
        self.skip_block(3)
        # stator 1
        self.skip_block(2)
# ---
        self.skip_block(62)

        ANZ_FORCE_AREAS = self.next_block("i")[0]

        if ANZ_FORCE_AREAS > 3:
            ANZ_FORCE_AREAS = 3

        self.skip_block()
        self.skip_block(2 * ANZ_FORCE_AREAS)
        #        self.skip_block(14)
        self.skip_block(10)
        self.delta_node_angle = self.next_block("f")[1]  # rad
        self.skip_block(3)
        self.skip_block(2 * 3 + 6 * 100 * 3)
        self.skip_block(30)
        self.skip_block(11 * 4)  # supel_EMODUL(I) ..
        self.skip_block()  # NODE_ISA_NODE_REC_ND_SV_PNTR
        self.skip_block(1 * 4)  # slot_insulation
        self.skip_block()  # NOM_CURRENT
        # PR_BASIC_LOSS_DATA ..
        # TOT_MAGNET_AREA
        # MOVE_EXTERN
        # MOVE_ARMATURE
        self.skip_block(4)
        try:
            self.pole_pairs, self.poles_sim = self.next_block("i")[:2]
        except:
            pass
        self.SLOT_WIRE_DIAMETER = self.next_block("f")
        self.SLOT_WIRE_NUMBERS = self.next_block("i")
        self.skip_block(20*(3 + 2 * 20))  # BASE_FREQUENCY ..
        self.skip_block(2)  # R_TORQUE .. NUM_NOLOAD_EX_CURRENT_STEPS
        try:
            (self.R_CURRENT,
             self.R_LOAD_VOLTAGE,
             self.R_NOLOAD_VOLTAGE) = self.next_block("f")
        except:
            pass
        x = self.next_block("f")
        self.R_COSPHI = x[0]
        self.R_BETA_OPT = x[1:]
        self.skip_block(10)  # R_FLUX_LOAD. NUM_NOLOAD_EX_CURRENT_STEPS

        if (FC_NUM_MOVE_LOSSES > 2 and NUM_FE_EVAL_MOVE_STEP > 1
                and FC_NUM_BETA_ID > 1):
            self.el_fe_induction_1.append([[]])
            self.el_fe_induction_2.append([[]])
            self.eddy_cu_vpot.append([[]])
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                b = self.next_block("h")
                logger.debug("el_fe_induction move losses (2) step %d: %d",
                            i, len(b))
                self.el_fe_induction_1[2][0].append(b)
                self.el_fe_induction_2[2][0].append(self.next_block("h"))
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.eddy_cu_vpot[2][0].append(self.next_block("h"))

        self.skip_block()
        self.skip_block(2 * 3)  # MAX_LOSS_EVAL_STEPS
        try:
            self.Q_SLOTS_NUMBER, self.M_PHASE_NUMBER = self.next_block("i")[:2]
        except:
            pass
        try:
            self.N_LAYERS_SLOT, self.N_WIRES_PER_SLOT = self.next_block("i")[
                :2]
        except:
            pass

        self.skip_block(1)
        self.skip_block(10 * 100)  # num_index_cad
        self.skip_block(1 * 100)
        self.skip_block()  # index_cad
        self.skip_block(1 * 4)  # heat_tranfer_coeff
        self.skip_block(2 * 2)
        self.skip_block()
        self.skip_block(2 * 4)
        self.skip_block(3)
        self.skip_block(1 * 64)  # bnodes_mech
        self.skip_block(5)
        self.SUPEL_ISA_SUPEL_REC_SE_FILLFACTOR = self.next_block("f")
        self.ELEM_ISA_ELEM_REC_LOSS_DENS = self.next_block("f")
        self.skip_block(3)
        self.skip_block(1 * 64)
        self.ROTOR_CUR_EXIST = self.next_block("?")[0]
        self.skip_block(20)  # mcmax = 20
        self.skip_block(4)
        self.NUM_SE_MAGN_KEYS = self.next_block("i")[0]

    def next_block(self, fmt):
        """
        Read binary data and return unpacked values according to format string.

        Arguments:
            fmt: Format string (see python struct module)
        """

        fmt_ = fmt.replace("?", "i")

        blockSize = struct.unpack_from("=i", self.file, self.pos)[0]
        self.pos += 4
        unpacked = []
        try:
            unpacked = struct.iter_unpack("=" + fmt_,
                                          self.file[self.pos:self.pos
                                                    + blockSize])
            unpacked = [x for x in unpacked]

        except AttributeError:  # python 2 has no iter_unpack
            chunksize = struct.calcsize("=" + fmt_)
            offset = self.pos
            for j in range(blockSize // chunksize):
                unpacked.append(struct.unpack_from("=" + fmt_,
                                                   self.file,
                                                   offset))
                offset += chunksize
            logger.debug("%s: %d %d", fmt_, blockSize, len(unpacked))
        except struct.error as e:
            logger.warning("Invalid Blocksize %s at pos %i",
                           blockSize, self.pos-4)

        self.pos += blockSize + 4

        fmt_ = ""
        for s in re.findall(r"[0-9]*.|[0-9]*\?", fmt):
            if len(s) > 1 and s[-1] != "s":
                fmt_ += int(s[:-1]) * s[-1]
            else:
                fmt_ += s
        values = []
        for i, dtype in enumerate(re.findall(r"\?|[0-9]*s?", fmt_)[:-1]):
            if dtype == "?":
                values.append([bool(u[i]) for u in unpacked])
            elif "s" in dtype:
                values.append([u[i].decode('latin-1') for u in unpacked])
            else:
                values.append([u[i] for u in unpacked])

        if len(fmt) == 1:
            return values[0]
        else:
            return values

    def skip_block(self, skips=1):
        """
        Proceed to the next block without reading any data.

        Arguments:
            skips: number of blocks to be skipped
        """
        while skips > 0:
            blockSize = struct.unpack_from("=i", self.file, self.pos)[0]
            self.pos += 4 + blockSize + 4
            skips -= 1


class Isa7(object):
    """
    The ISA7 Femag model

    Arguments:
        filename: name of I7/ISA7 file to be read
    """

    color = {1: [1.0, 0.0, 0.0],  # RED
             2: [0.0, 1.0, 0.0],  # GREEN
             3: [1.0, 0.8, 0.0],  # DARKYELLOW
             4: [0.0, 0.5019607843137255, 1.0], # BLUE
             5: [0.9803921568627451, 0.0, 1.0], # MAGENTA
             6: [0.0, 1.0, 0.8235294117647058], # CYAN
             7: [1.0, 1.0, 1.0],  # BRIGHTWHITE
             8: [0.0, 0.0, 0.0],  # BLACK
             9: [0.0, 0.0, 0.5882352941176471], # DARKBLUE
             10: [0.6666666666666666, 0.0, 0.0],# DARKRED
             11: [0.6666666666666666, 1.0, 0.0],# DARKGREEN
             12: [1.0, 0.6274509803921569, 0.0],# BROWN
             13: [0.0, 0.0, 1.0], # MEDIUMBLUE
             14: [0.6666666666666666, 0.0, 1.0], # VIOLET
             15: [0.0, 0.8235294117647058, 1.0], # SKYBLUE
             16: [0.8274509803921568, 0.8274509803921568, 0.8274509803921568]} # LIGHGREY

    def __init__(self, reader):
        try:
            self.state_of_problem = reader.state_of_problem
        except:
            pass
        self.points = [Point(x, y)
                       for x, y in zip(reader.POINT_ISA_POINT_REC_PT_CO_X,
                                       reader.POINT_ISA_POINT_REC_PT_CO_Y)]

        self.lines = [Line(self.points[abs(pk1) - 1], self.points[abs(pk2) - 1])
                      for pk1, pk2 in zip(reader.LINE_ISA_LINE_REC_LN_PNT_1,
                                          reader.LINE_ISA_LINE_REC_LN_PNT_2)]
        logger.debug("Nodes")
        self.nodes = [
            Node(n + 1,
                 reader.NODE_ISA_NODE_REC_ND_BND_CND[n],
                 reader.NODE_ISA_NODE_REC_ND_PER_NOD[n],
                 reader.NODE_ISA_ND_CO_RAD[n],
                 reader.NODE_ISA_ND_CO_PHI[n],
                 reader.NODE_ISA_NODE_REC_ND_CO_1[n],
                 reader.NODE_ISA_NODE_REC_ND_CO_2[n],
                 reader.NODE_ISA_NODE_REC_ND_VP_RE[n],
                 reader.NODE_ISA_NODE_REC_ND_VP_IM[n])
            for n in range(len(reader.NODE_ISA_NODE_REC_ND_BND_CND))]

        logger.debug("Nodechains")
        self.nodechains = []
        for nc in range(len(reader.NDCHN_ISA_NDCHN_REC_NC_NOD_1)):
            nd1 = reader.NDCHN_ISA_NDCHN_REC_NC_NOD_1[nc]
            nd2 = reader.NDCHN_ISA_NDCHN_REC_NC_NOD_2[nc]
            ndm = reader.NDCHN_ISA_NDCHN_REC_NC_NOD_MID[nc]
            try:
                node1 = self.nodes[abs(nd1) - 1]
                nodem = self.nodes[ndm - 1]
                node2 = self.nodes[abs(nd2) - 1]

                if nd1 < 0 or nd2 < 0:
                    nodes = node1, nodem, node2
                elif ndm > 0:
                    nodes = node1, nodem, node2
                else:
                    nodes = node1, None, node2

                self.nodechains.append(
                    NodeChain(nc + 1, nodes))
            except IndexError:
                logger.warning('IndexError in nodes')
                raise  # preserve the stack trace

        self.elements = []
        logger.debug("Elements")
        for e in range(len(reader.ELEM_ISA_EL_NOD_PNTR)):
            ndkeys = []
            ndk = reader.ELEM_ISA_EL_NOD_PNTR[e]

            while ndk > 0:
                ndkeys.append(reader.ELE_NOD_ISA_ND_KEY[ndk - 1])
                ndk = reader.ELE_NOD_ISA_NXT_ND_PNTR[ndk - 1]

            vertices = [self.nodes[k - 1] for k in ndkeys]
            try:
                loss_dens = reader.ELEM_ISA_ELEM_REC_LOSS_DENS[e]
            except (IndexError, AttributeError):
                loss_dens = 0
            try:
                temperature = reader.ELEM_ISA_ELEM_REC_TEMPERATURE[e]
            except (IndexError, AttributeError):
                temperature = 20
            self.elements.append(
                Element(e + 1,
                        ElType(reader.ELEM_ISA_ELEM_REC_EL_TYP[e]),
                        reader.ELEM_ISA_ELEM_REC_EL_SE_KEY[e] - 1,
                        vertices,
                        (reader.ELEM_ISA_ELEM_REC_EL_RELUC[e],
                         reader.ELEM_ISA_ELEM_REC_EL_RELUC_2[e]),
                        (reader.ELEM_ISA_ELEM_REC_EL_MAG_1[e],
                         reader.ELEM_ISA_ELEM_REC_EL_MAG_2[e]),
                        loss_dens,  # in W/m³
                        reader.BR_TEMP_COEF/100,
                        temperature)   # in 1/K
            )

        logger.debug("SuperElements")
        self.superelements = []
        for se in range(len(reader.SUPEL_ISA_SE_NDCHN_PNTR)):
            nc_keys = []
            nc_ptr = reader.SUPEL_ISA_SE_NDCHN_PNTR[se]

            while nc_ptr > 0:
                nc_keys.append(reader.SE_NDCHN_ISA_NC_KEY[nc_ptr - 1])
                nc_ptr = reader.SE_NDCHN_ISA_NXT_NC_PNTR[nc_ptr - 1]

            nodechains = []
            for nck in nc_keys:
                if nck > 0:
                    nodechains.append(self.nodechains[abs(nck) - 1])
                else:
                    nodechains.append(self.nodechains[abs(nck) - 1].reverse())

            el_keys = []
            el_ptr = reader.SUPEL_ISA_SE_EL_PNTR[se]

            while el_ptr > 0:
                el_keys.append(reader.SE_EL_ISA_EL_KEY[el_ptr - 1])
                el_ptr = reader.SE_EL_ISA_NXT_EL_PNTR[el_ptr - 1]

            elements = [self.elements[elk - 1]
                        for elk in el_keys]
            try:
                fillfactor = reader.SUPEL_ISA_SUPEL_REC_SE_FILLFACTOR[se]
            except:
                fillfactor = 1
            try:
                temp_coef = reader.SUPEL_ISA_SUPEL_REC_SE_TEMP_COEF[se]
            except:
                temp_coef = 0
            try:
                temperature = reader.SUPEL_ISA_SUPEL_REC_SE_TEMPERATURE[se]
            except:
                temperature = 20
            self.superelements.append(
                SuperElement(se + 1,
                             reader.SUPEL_ISA_SUPEL_REC_SE_SR_KEY[se] - 1,
                             elements,
                             nodechains,
                             reader.SUPEL_ISA_SUPEL_REC_SE_COL[se],
                             nc_keys,
                             reader.SUPEL_ISA_SUPEL_REC_SE_MCV_TYP[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_COND_TYP[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_CONDUC[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_LENGHT[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_VEL_SYS[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_VELO_1[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_VELO_2[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_CURD_RE[se],
                             reader.SUPEL_ISA_SUPEL_REC_SE_CURD_IM[se],
                             fillfactor, temp_coef, temperature))

        logger.debug("Subregions")
        self.subregions = []
        for sr in range(len(reader.SR_ISA_SR_SE_PNTR)):
            se_keys = []
            se_ptr = reader.SR_ISA_SR_SE_PNTR[sr]

            while se_ptr > 0:
                se_keys.append(reader.SR_SE_ISA_SE_KEY[se_ptr - 1])
                se_ptr = reader.SR_SE_ISA_NXT_SE_PNTR[se_ptr - 1]

            superelements = [self.superelements[sek - 1]
                             for sek in se_keys]

            nodechains = []
            nc_keys = []
            for se in superelements:
                nc_keys.extend([abs(nc.key) for nc in se.nodechains])
            nc_keys = [nck for nck, count
                       in Counter(nc_keys).items() if count < 2]
            for se in superelements:
                nodechains.extend([nc
                                   for nc in se.nodechains
                                   if abs(nc.key) in nc_keys])

            self.subregions.append(
                SubRegion(sr + 1,
                          reader.SR_ISA_SR_REC_SR_TYP[sr],
                          reader.SR_ISA_SR_REC_SR_COL[sr],
                          reader.SR_ISA_SR_REC_SR_NAME[sr].strip(),
                          reader.SR_ISA_SR_REC_SR_NTURNS[sr],
                          reader.SR_ISA_SR_REC_SR_CUR_DIR[sr],
                          reader.SR_ISA_SR_REC_SR_WB_KEY[sr] - 1,
                          superelements,
                          nodechains))

        logger.debug("Windings")
        self.windings = []
        try:
            for wd in range(len(reader.WB_ISA_WB_SR_PNTR)):
                sr_keys = []
                sr_ptr = reader.WB_ISA_WB_SR_PNTR[wd]

                while sr_ptr > 0:
                    sr_keys.append(reader.WB_SR_ISA_SR_KEY[sr_ptr - 1])
                    sr_ptr = reader.WB_SR_ISA_NXT_SR_PNTR[sr_ptr - 1]

                subregions = [self.subregions[srk - 1]
                              for srk in sr_keys]

                self.windings.append(
                    Winding(wd + 1,
                            reader.WB_ISA_WB_REC_WB_NAME[wd],
                            subregions,
                            reader.WB_ISA_WB_REC_WB_TURN[wd],
                            reader.WB_ISA_WB_REC_WB_GCUR_RE[wd],
                            reader.WB_ISA_WB_REC_WB_GCUR_IM[wd],
                            reader.WB_ISA_WB_REC_WB_IMPDZ_RE[wd],
                            reader.WB_ISA_WB_REC_WB_IMPDZ_IM[wd],
                            reader.WB_ISA_WB_REC_WB_VOLT_RE[wd],
                            reader.WB_ISA_WB_REC_WB_VOLT_IM[wd]))
        except:
            pass
        logger.info("Total nodes %d elements %d superelements %d subregions %d",
                    len(self.nodes), len(self.elements),
                    len(self.superelements),
                    len(self.subregions))

        # positions of all elements
        self.element_pos = np.array([e.center
                                     for e in self.elements])

        for a in ('FC_RADIUS', 'pole_pairs', 'poles_sim', 'move_action',
                  'layers', 'coil_span', 'delta_node_angle', 'speed',
                  'MAGN_TEMPERATURE', 'BR_TEMP_COEF',
                  'MA_SPEZ_WEIGHT', 'CU_SPEZ_WEIGHT'):
            try:
                setattr(self, a, getattr(reader, a))
            except AttributeError:
                pass
        if getattr(reader, 'pole_pairs', 0):
            self.num_poles = 2*self.pole_pairs
        if getattr(reader, 'slots', 0):
            self.num_slots = reader.slots
        try:
            self.arm_length = reader.arm_length*1e-3  # in m
        except AttributeError:
            # missing arm_length
            pass

        self.airgap_inner_elements = [] # used for rotate
        self.airgap_outer_elements = []
        if getattr(self, 'FC_RADIUS', 0) > 0:  # Note: cosys r/phi only
            # TODO: handle multiple airgaps
            airgap_center_elements = []
            for n in self.nodes:
                n.outside = np.linalg.norm(n.xy) > self.FC_RADIUS

            for e in self.elements:
                outside = [v.outside
                           for v in e.vertices]
                if all(outside):
                    self.airgap_outer_elements.append(e)
                elif any(outside):
                    airgap_center_elements.append(e)
                else:
                    self.airgap_inner_elements.append(e)

            self.airgap_center_elements = sorted(airgap_center_elements,
                                                 key=lambda e: np.arctan2(
                                                     e.center[1], e.center[0]))
        else:  # assume a linear machine
            # TODO read and check pole_width
            airgap_positions = []
            self.airgap_center_elements = []
            nxy = np.array([n.xy for n in self.nodes])
            # width and height of model:
            w, h = np.max(nxy[:, 0])-np.min(nxy[:, 0]), np.max(nxy[:, 1])-np.min(nxy[:, 1])
            for se in self.superelements:
                if (se.subregion is None and se.elements and not (
                        se.elements[0].is_magnet() or
                        se.elements[0].is_lamination()) and
                    se.elements[0].el_type == ElType.LinearRectangle):
                    axy = np.array([n.xy for e in se.elements
                                    for n in e.vertices])
                    # width and height of superelement
                    sew = np.max(axy[:, 0])-np.min(axy[:, 0])
                    seh = np.max(axy[:, 1])-np.min(axy[:, 1])
                    if np.isclose(sew, w) or np.isclose(seh, h):
                        horiz = sew > seh
                        if horiz:
                            airgap_positions.append(axy[0][1])
                        else:
                            airgap_positions.append(axy[0][0])
                        self.airgap_center_elements.append(se.elements)

            if self.airgap_center_elements:
                # TODO check airgap center
                self.airgap_center_elements = self.airgap_center_elements[-1]
                if horiz:
                    airgap_positions.append(np.min(nxy[:, 1]))
                else:
                    airgap_positions.append(np.max(nxy[:, 0]))
            if airgap_positions:
                amin, amax = min(airgap_positions), max(airgap_positions)

                def is_outside(n):
                    if horiz:
                        return (n.y > amax or n.y < amin)
                    return (n.x > amax or n.x < amin)

                for n in self.nodes:
                    n.outside = is_outside(n)

                for e in self.elements:
                    outside = [v.outside
                               for v in e.vertices]
                    if all(outside):
                        self.airgap_outer_elements.append(e)
                    elif not all(outside):
                        self.airgap_inner_elements.append(e)

        for se in self.superelements:
            if se.nodechains:
                se.outside = se.nodechains[0].node1.outside

        self.pos_el_fe_induction = np.asarray(reader.pos_el_fe_induction)
        if self.pos_el_fe_induction.shape[0]:
            # pos_el_fe_induction: strictly monotone increasing sequence
            a = self.pos_el_fe_induction
            a = np.concatenate(([a[0]], a[1:][a[1:] > a[:-1]]))
            if a.shape[0] < self.pos_el_fe_induction.shape[0]:
                self.pos_el_fe_induction = a
        try:
            self.beta_loss = np.asarray(reader.beta_loss)
            self.curr_loss = np.array([c/np.sqrt(2) for c in reader.curr_loss])
            logger.debug("Beta loss %s curr loss %s",
                         reader.beta_loss, reader.curr_loss)
        except AttributeError:
            try:
                self.iq = np.asarray(reader.iq)
                self.id = np.array(reader.id)
                logger.debug("Id loss %s iq %s",
                             reader.id, reader.iq)
            except AttributeError:
                pass

        try:
            self.current_id = np.asarray(reader.CURRENT_ID)
            self.ide_flux = np.asarray(reader.IDE_FLUX)
            self.ide_beta = np.asarray(reader.IDE_BETA)
        except AttributeError:
            pass

        try:
            flx_fac = 1000
            if isinstance(reader.el_fe_induction_1, list):
                pass
            else:
                if reader.el_fe_induction_1.dtype == 'int16':
                    flx_fac = 1000
                else:
                    flx_fac = 1
            el_fe_ind = [np.array(reader.el_fe_induction_1).T/flx_fac,
                        np.array(reader.el_fe_induction_2).T/flx_fac]
            eddy_cu_vpot = np.array(reader.eddy_cu_vpot).T/1000
            if len(el_fe_ind[0].shape) == 4:
                pdim = self.pos_el_fe_induction.shape[0]
                if pdim < el_fe_ind[0].shape[1]:
                    el_fe_ind = [el_fe_ind[0][:,:pdim, :, :],
                                 el_fe_ind[1][:,:pdim, :, :]]
                    eddy_cu_vpot = eddy_cu_vpot[:,:pdim, :, :]
        except (ValueError, TypeError) as e:
            # inhomogenous array
            l = len(reader.el_fe_induction_1[0][0])
            shape = []
            for i in reader.el_fe_induction_1:
                for j in i:
                    n = 0
                    for k in j:
                        if len(k) < l:
                            break
                        n += 1
                    if n > 0:
                        shape.append(n)

            el_fe_ind = [np.array([[reader.el_fe_induction_1[0][0][:shape[0]]]]).T/flx_fac,
                        np.array([[reader.el_fe_induction_2[0][0][:shape[0]]]]).T/flx_fac]
            eddy_cu_vpot = np.array([[reader.eddy_cu_vpot[0][0][:shape[0]]]]).T/1000

        self.el_fe_induction_1 = el_fe_ind[0]
        self.el_fe_induction_2 = el_fe_ind[1]
        self.eddy_cu_vpot = eddy_cu_vpot

        self.PS_FILFACTOR_CU = reader.PS_FILFACTOR_CU
        self.PS_LENGTH_CU = reader.PS_LENGTH_CU
        self.PS_SIGMA_CU = reader.PS_SIGMA_CU

        self.iron_loss_coefficients = getattr(
            reader, 'iron_loss_coefficients', [])

    def get_subregion(self, name):
        """return subregion by name"""
        for s in self.subregions:
            if s.name == name.strip():
                return s
        raise ValueError('no such subregion "{}" in this model'.format(name))

    def wdg_elements(self):
        """return elements in winding region"""
        return [el for el in self.elements
                if el.superelement.condtype != 0]

    def magnet_super_elements(self):
        """return superelements which are magnets"""
        return [self.superelements[i]
                for i in set([el.se_key for el in self.magnet_elements()])]

    def magnet_elements(self):
        """return elements which are magnets"""
        return [e for e in self.elements if e.is_magnet()]

    def get_element(self, x, y):
        """return element at pos x,y"""
        k = np.argmin(np.linalg.norm(self.element_pos - (x, y), axis=1))
        return self.elements[k]

    def get_super_element(self, x, y):
        """return superelement at pos x,y"""
        e = self.get_element(x, y)
        try:
            return [s for s in self.superelements
                    if e.key in [se.key for se in s.elements]][0]
        except IndexError:
            return None

    def flux_density(self, el, icur, ibeta):
        """return move pos and flux density in model coordinates (bx, by)
        of element for current and beta

        Arguments:
          el: element
          icur, ibeta: current and beta (load) index (0: noload, 1: zero, 2: as specified)
        """
        ekey = el.key-1
        b1 = np.array(self.el_fe_induction_1[ekey, :, icur, ibeta])
        b2 = np.array(self.el_fe_induction_2[ekey, :, icur, ibeta])
        return dict(
            pos=self.pos_el_fe_induction,
            bx=b1,
            by=b2)

    def get_areas(self):
        """ return areas (in m²) of inner and outer regions slots, iron, magnets
        """
        try:
            return self.areas
        except AttributeError:
            self.areas = [{'iron': 0, 'slots': 0, 'magnets': 0},
                          {'iron': 0, 'slots': 0, 'magnets': 0}]
        scf = self.scale_factor()
        for se in self.superelements:
            r = 0 if se.outside else 1
            if se.subregion:
                if se.subregion.winding:
                    self.areas[r]['slots'] += scf*se.area()
                    continue

            if se.mcvtype or se.elements[0].is_lamination():
                self.areas[r]['iron'] += se.area()*scf
            else:
                a = [e.area for e in se.elements if e.is_magnet()]
                self.areas[r]['magnets'] += sum(a)*scf
        return self.areas

    def get_mass(self):
        """ return mass (in kg) of material conductors, iron, magnets
        """
        try:
            return self.mass
        except AttributeError:
            self.mass = [{'iron': 0, 'conductors': 0, 'magnets': 0},
                         {'iron': 0, 'conductors': 0, 'magnets': 0}]
        scf = self.scale_factor()
        for se in self.superelements:
            r = 0 if se.outside else 1
            if se.subregion:
                if se.subregion.winding:
                    spw = self.CU_SPEZ_WEIGHT*1e3
                    l = self.PS_LENGTH_CU[r]*1e-2
                    slf = self.PS_FILFACTOR_CU[r]
                    self.mass[r]['conductors'] += scf*se.area()*self.arm_length*spw*l*slf
                    continue

            if se.mcvtype or se.elements[0].is_lamination():
                try:
                    spw = self.iron_loss_coefficients[se.mcvtype-1][
                        'spec_weight']*1e3  # kg/m³
                    fillfact = self.iron_loss_coefficients[se.mcvtype-1][
                        'fillfactor']
                except IndexError:
                    spw = 7.8e3
                    fillfact = 1
                    #logger.warning('missing iron losscoeffs using defaults')
                m = scf*self.arm_length*se.area()*spw*fillfact
                self.mass[r]['iron'] += m
            else:
                try:
                    spw = self.MA_SPEZ_WEIGHT*1e3
                    a = [e.area for e in se.elements if e.is_magnet()]
                    self.mass[r]['magnets'] += sum(a)*scf*self.arm_length*spw
                except AttributeError:
                    pass
        return self.mass

    def get_iron_subregions(self) -> list:
        """return names of subregions with lamination

        Returns:
          subregs: list of str

        """
        return [sr.name for sr in self.subregions
                if sr.superelements[0].mcvtype or sr.superelements[0].elements[0].is_lamination()]

    def _axis_ratio(self, apos, br, bt):
        from .utils import fft
        pulsating = 0
        brtmax = np.max(br-np.mean(br)), np.max(bt-np.mean(bt))
        if np.all(np.isclose(brtmax, 0)) or np.any(np.isclose(brtmax, 0)):
            return pulsating

        br0 = fft(apos, br-np.mean(br))
        x = br0['a']*np.cos(2*np.pi*apos/br0['T0']+br0['alfa0'])
        bt0 = fft(apos, bt-np.mean(br))
        y = bt0['a']*np.cos(2*np.pi*apos/bt0['T0']+bt0['alfa0'])
        if (br0['a'] > br0['nue'][self.pole_pairs]
            or bt0['a'] > bt0['nue'][self.pole_pairs]):
            return pulsating

        kmax = np.argmax(np.linalg.norm((x, y), axis=0))
        kmin = np.argmin(np.linalg.norm((x, y), axis=0))
        a = np.linalg.norm((x[kmax], y[kmax]))
        b = np.linalg.norm((x[kmin], y[kmin]))
        return b/a  # ecc: np.sqrt(1-b**2/a**2))


    def calc_iron_loss(self, icur: int, ibeta: int, pfefun, bmin=0.1) -> dict:
        """ calculate iron loss using last simulation results

        Args:
          icur: current index
          ibeta: beta index
          pfefun: custom function with parameters:
            Bxnu, Bynu (1d array): flux density values in T
            fnu (1d array): frequency values in Hz
            losscoeffs (dict): material properties
            axr: float (optional)
          bmin (float): lower limit of flux density amplitudes

        Returns:
          loss values name of subregion (string),
            hysteresis loss, eddy current loss, exc loss

        """
        from .utils import fft
        from inspect import signature
        # check if axis ratio is needed
        need_axratio = len(signature(pfefun).parameters) > 4
        # eliminate double values at end
        # TODO: adapt to linear machines

        pos = [p
               for p in self.pos_el_fe_induction
               if p < 2*np.pi/self.pole_pairs] + [2*np.pi/self.pole_pairs]
        apos = np.array(pos)/np.pi*180
        i = len(pos)
        sreg = {}
        f1 = self.speed/60
        scf = self.scale_factor()
        for sr in [self.get_subregion(sname)
                   for sname in self.get_iron_subregions()]:
            losses = []
            for se in sr.superelements:
                spw = self.iron_loss_coefficients[se.mcvtype-1]['spec_weight']
                fillfact = self.iron_loss_coefficients[se.mcvtype-1]['fillfactor']
                for e in se.elements:
                    br = self.el_fe_induction_1[e.key-1, 0:i+1, icur, ibeta]
                    bt = self.el_fe_induction_2[e.key-1, 0:i+1, icur, ibeta]
                    b1 = fft(apos, br, pmod=2)
                    b2 = fft(apos, bt, pmod=2)
                    if need_axratio:
                        axr = self._axis_ratio(apos, br, bt)
                    blen = max(len(b1['nue']), len(b2['nue']))
                    if len(b1['nue']) < blen:
                        b1['nue'] = b1['nue'] + [0]*(blen-len(b1['nue']))
                    if len(b2) < blen:
                        b2['nue'] = b2['nue'] + [0]*(blen-len(b2['nue']))
                    bnxy = np.array(
                        [(n, b[0], b[1]) for n, b in enumerate(
                            zip(b1['nue'], b2['nue']))
                         if b[0] > bmin or b[1] > bmin]).T
                    if bnxy.size > 0:
                        fnu = np.array([f1*nue for nue in bnxy[0]])
                        if need_axratio:
                            phy, pec, pex = pfefun(
                                bnxy[1]/fillfact, bnxy[2]/fillfact, fnu,
                            self.iron_loss_coefficients[se.mcvtype-1], axr)
                        else:
                            phy, pec, pex = pfefun(
                                bnxy[1]/fillfact, bnxy[2]/fillfact, fnu,
                                self.iron_loss_coefficients[se.mcvtype-1])

                        pl = [1e3*spw*e.area*l
                              for l in (sum(phy), sum(pec), sum(pex))]
                        losses.append(pl)
                    else:
                        logger.debug("Empty %s, %s", b1, b2)
            logger.debug("%s: %s", sr.name, losses)
            if losses:
                sreg[sr.name] = (scf*self.arm_length*np.sum(losses, axis=0)).tolist()
            else:
                sreg[sr.name] = [0,0,0]
        return sreg

    def get_minmax_temp(self):
        def node_subregion(subregion_name):
            node_temperature = []
            for i in self.subregions:
                if i.name.lower() == subregion_name:
                    for j in i.elements():
                        for k in j.vertices:
                            node_temperature.append(k.vpot[-1])
            ndtemp = np.unique(node_temperature)
            return [np.amax(ndtemp), np.amin(ndtemp), np.mean(ndtemp)]

        zero_temp = [0.0, 0.0, 0.0]
        component_temperature = dict(styoke=zero_temp, stteeth=zero_temp,
                                     magnet=zero_temp, winding=zero_temp)
        for i in self.subregions:
            if i.name:
                sreg_name = i.name.lower()
                if sreg_name in ('stza', 'stth', 'stteeth'):
                    component_temperature['stteeth'] = node_subregion(sreg_name)
                elif sreg_name in ('styk', 'styoke', 'stjo'):
                    component_temperature['styoke'] = node_subregion(sreg_name)
                elif sreg_name == 'iron':
                    component_temperature['iron'] = node_subregion(sreg_name)
                elif sreg_name == 'pmag':
                    component_temperature['magnet'] = node_subregion(sreg_name)
                else:
                    pass

        wdg_temps = []
        for j in self.wdg_elements():
            for k in j.vertices:
                wdg_temps.append(k.vpot[-1])
        wdg_temp = np.unique(wdg_temps)
        component_temperature['winding'] = [np.amax(wdg_temp),
                                            np.amin(wdg_temp),
                                            np.mean(wdg_temp)]
        return component_temperature

    def flux_dens(self, x, y, icur, ibeta):
        el = self.get_element(x, y)
        return self.flux_density(el, icur, ibeta)

    def demagnetization(self, el, icur, ibeta):
        """return demagnetization Hx, Hy at element
        Arguments:
          el: element
          icur, ibeta: current, beta index"""
        flxdens = self.flux_density(el, icur, ibeta)
        return (flxdens['pos'], el.demag_b((flxdens['bx'], flxdens['by']),
                                           self.MAGN_TEMPERATURE))

    def demag_situation(self, icur, ibeta, hlim):
        """return h max, h avg, area, pos for demag situation for
        each magnet
        Arguments:
          icur: cur amplitude index
          ibeta: beta angle index (load)
          hlim: limit of demagnetization (kA/m)
        """
        results = []
        for se in self.magnet_super_elements():
            elements = np.array(se.elements)
            demag = np.array([self.demagnetization(el, icur, ibeta)[1]
                              for el in elements])
            ind = np.unravel_index(np.argmax(demag, axis=None), demag.shape)
            dmax = demag[:, ind[1]]
            area_tot = np.sum([e.area for e in elements])
            area_demag = np.sum([e.area for e in elements[-dmax < hlim]])
            results.append(dict(
                h_max=-demag[ind],
                h_avg=-np.average(dmax),
                area_tot=area_tot,
                area_demag=area_demag,
                pos=self.pos_el_fe_induction[ind[1]]))
        return results

    def rotate(self, alpha):
        if alpha:
            for n in self.nodes:
                if not n.outside:
                    n.xy = (np.cos(alpha)*n.x -np.sin(alpha)*n.y,
                            np.sin(alpha)*n.x + np.cos(alpha)*n.y)
            self.elements = self.airgap_outer_elements + self.airgap_inner_elements
        else: # reset rotation
            for n in self.nodes:
                if not n.outside:
                    n.xy = n.x, n.y
            self.elements = (self.airgap_outer_elements + self.airgap_center_elements +
                             self.airgap_inner_elements)

    def scale_factor(self):
        '''Returns the scale factor
        Parameters
        ----------
        None

        Returns
        -------
        scale_factor : int
        '''
        try:
            poles = 2*self.pole_pairs
        except:
            poles = 2
        try:
            poles_sim = self.poles_sim
        except:
            poles_sim = poles

        scale_factor = poles/poles_sim
        return scale_factor

    def get_magnet_data(self, ibeta=0, icur=0) -> list:
        '''Extract magnet data from nc file

        Args:
            nc: nc object
            icur, ibeta: load case

        Returns:
          pm_data: list of magnet data
        '''
        mag_spels = self.magnet_super_elements()
        if len(mag_spels) / self.poles_sim == 1:
            mag_spels = [mag_spels[0]]

        # prepare data for ialh method
        # conductivity and permeability of the magnets
        cond = 0
        mur = 0
        # read boundary nodes
        for se in mag_spels:
            cond = se.conduc
            if cond <= 0:
                cond = 625000
                logging.warning('Magnet conductivity  <= 0, using 625000 S/m')
            mur = np.abs(1/se.elements[0].reluc[0])
            logging.debug('Magnet: mur=%s, conductivity=%s', mur, cond)

        # stationary case, no rotation
        poles = 0
        try:
            poles = self.num_poles
        except AttributeError:
            pass

        if poles == 0:  # no rotation
            freq = self.speed
            time_vec = np.linspace(0, 1/freq, len(self.pos_el_fe_induction))
            pos = dict(time=time_vec.tolist(),
                       freq=freq,
                       t=float(1/freq))
            # reset num.poles
            poles = 1
        else:
            if self.move_action == 0:
                phi = self.pos_el_fe_induction*180/np.pi
                pos = dict(phi=phi,
                           speed=self.speed)
            else:
                pos = dict(displ=self.pos_el_fe_induction,
                           speed=self.speed)
        # prep dictionary for the loss calculation
        pm_data = []
        for i, se in enumerate(mag_spels):
            ecp = [e.center for e in se.elements]
            geometry = se.get_rect_geom()

            bxy = []
            for e in se.elements:
                theta = np.arctan2(float(e.center[1]),
                                   float(e.center[0]))
                fd = self.flux_density(e, icur, ibeta)
                bxy.append(Trot(-theta).dot((fd['bx'], fd['by'])))
            #= np.moveaxis(bxy, 1, 0)
            pd = dict(name='pm_data_se' + str(se.key),
                      hm=geometry['h'],
                      wm=geometry['w'],
                      lm=self.arm_length,
                      alpha=geometry['alpha'],
                      ls=self.arm_length,
                      sigma=cond,
                      mur=mur,
                      loadcase=ibeta,
                      numpoles=poles,
                      bl=transform_flux_density(geometry['alpha'],
                                                np.array(bxy)),
                      elcp=transform_coord(geometry, ecp),
                      area=se.area(),
                      spel_key=se.key)
            pd.update(pos)

            pm_data.append(pd)
        return pm_data


class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.xy = x, y


class Line(object):
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2


class BaseEntity(object):
    def __init__(self, key):
        self.key = key


class Node(BaseEntity):
    def __init__(self, key, bndcnd, pernod, r, phi, x, y, vpot_re, vpot_im):
        super(self.__class__, self).__init__(key)
        self.bndcnd = bndcnd
        self.pernod = pernod
        self.r = r
        self.phi = phi
        self.x = x
        self.y = y
        self.xy = x, y
        self.vpot = vpot_re, vpot_im
        self.outside = True

    def on_boundary(self):
        return self.bndcnd != 0 or self.pernod != 0


class NodeChain(BaseEntity):
    def __init__(self, key, nodes):
        super(self.__class__, self).__init__(key)
        self.node1 = nodes[0]
        self.nodemid = nodes[1]
        self.node2 = nodes[2]
        if nodes[1] is None:
            self.nodes = (nodes[0], nodes[2])
        else:
            self.nodes = (nodes[0], nodes[1], nodes[2])

    def reverse(self):
        return NodeChain(self.key * (-1),
                         [self.node2, self.nodemid, self.node1])


class Element(BaseEntity):
    def __init__(self, key, el_type,
                 se_key, vertices, reluc, mag, loss_density,
                 br_temp_coef=0, temperature=20):
        super(self.__class__, self).__init__(key)
        self.el_type = el_type
        self.se_key = se_key
        self.vertices = vertices
        self.reluc = reluc
        self.mag = mag
        self.br_temp_coef = br_temp_coef
        self.loss_density = loss_density
        self.temperature = temperature
        if el_type == ElType.LinearTriangle:
            self.area = ((vertices[2].x - vertices[1].x) *
                         (vertices[0].y - vertices[1].y) -
                         (vertices[2].y - vertices[1].y) *
                         (vertices[0].x - vertices[1].x))/2
        elif el_type == ElType.LinearRectangle:
            self.area = ((vertices[2].x - vertices[1].x) *
                         (vertices[0].y - vertices[1].y) -
                         (vertices[2].y - vertices[1].y) *
                         (vertices[0].x - vertices[1].x) +
                         (vertices[3].x - vertices[2].x) *
                         (vertices[0].y - vertices[2].y) -
                         (vertices[3].y - vertices[2].y) *
                         (vertices[0].x - vertices[2].x))/2
        elif el_type == ElType.SquareTriangle:
            self.area = ((vertices[4].x - vertices[2].x) *
                         (vertices[0].y - vertices[2].y) -
                         (vertices[4].y - vertices[1].y) *
                         (vertices[0].x - vertices[2].x))/2
        elif el_type == ElType.SquareRectangle:
            self.area = ((vertices[4].x - vertices[2].x) *
                         (vertices[0].y - vertices[2].y) -
                         (vertices[4].y - vertices[2].y) *
                         (vertices[0].x - vertices[2].x) +
                         (vertices[6].x - vertices[4].x) *
                         (vertices[0].y - vertices[4].y) -
                         (vertices[6].y - vertices[4].y) *
                         (vertices[0].x - vertices[4].x))/2
        self.center = np.sum(
            [v.xy for v in vertices], axis=0)/len(vertices)

    def flux_density(self, cosys='cartes'):
        """return flux density components of this element converted to cosys: cartes, cylind, polar"""
        ev = self.vertices
        b1, b2 = 0, 0
        if self.el_type == ElType.LinearTriangle:
            y31 = ev[2].y - ev[0].y
            y21 = ev[1].y - ev[0].y
            x13 = ev[0].x - ev[2].x
            x21 = ev[1].x - ev[0].x
            a21 = ev[1].vpot[0] - ev[0].vpot[0]
            a31 = ev[2].vpot[0] - ev[0].vpot[0]
            delta = self.superelement.length * (y31 * x21 + y21 * x13)

            b1, b2 = ((x13 * a21 + x21 * a31) / delta,
                      (-y31 * a21 + y21 * a31) / delta)

        elif self.el_type == ElType.LinearRectangle:
            y31 = ev[2].y - ev[0].y
            y21 = ev[1].y - ev[0].y
            x13 = ev[0].x - ev[2].x
            x21 = ev[1].x - ev[0].x
            a21 = ev[1].vpot[0] - ev[0].vpot[0]
            a31 = ev[2].vpot[0] - ev[0].vpot[0]
            delta = self.superelement.length * (y31 * x21 + y21 * x13)
            b1_a = (x13 * a21 + x21 * a31) / delta
            b2_a = (y21 * a31 - y31 * a21) / delta

            y31 = ev[0].y - ev[2].y
            y21 = ev[3].y - ev[2].y
            x13 = ev[2].x - ev[0].x
            x21 = ev[3].x - ev[2].x
            a24 = ev[3].vpot[0] - ev[2].vpot[0]
            a34 = ev[0].vpot[0] - ev[2].vpot[0]
            delta = self.superelement.length * (y31 * x21 + y21 * x13)
            b1_b = (x13 * a24 + x21 * a34) / delta
            b2_b = (y21 * a34 - y31 * a24) / delta

            b1, b2 = ((b1_a + b1_b) / 2,
                      (b2_a + b2_b) / 2)

        if cosys == 'cartes':
            return (b1, b2)
        if cosys == 'polar':
            a = np.arctan2(self.center[1], self.center[0])
            br, bphi = np.array(((np.cos(a), np.sin(a)),
                                 (-np.sin(a), np.cos(a)))).dot(((b1), (b2)))
            return br, bphi
        if cosys == 'cylind':
            xm = np.sum([e.x for e in ev])
            rm = np.sum([e.vpot[0] for e in ev])
            if np.abs(xm) < 1e-6:
                rm = 0
            else:
                rm = rm/xm
            return -b1, -b2/rm

    def is_magnet(self):
        """return True if the element is a permanent magnet"""
        return abs(self.mag[0]) > 1e-5 or abs(self.mag[1]) > 1e-5

    def is_lamination(self):
        """return True if the element has lamination properties"""
        return self.reluc != (1.0, 1.0) and self.mag == (0.0, 0.0)

    def remanence(self, temperature=20):
        """return remanence Brx, Bry at element
        Arguments:
          temperature: (float) temperature in °C"""
        br_temp_corr = 1. + self.br_temp_coef*(temperature - 20.)
        return self.mag[0]*br_temp_corr, self.mag[1]*br_temp_corr

    def demagnetization(self, temperature=20):
        """return demagnetization Hx, Hy of this element"""
        return self.demag_b(self.flux_density(cosys='polar'), temperature)

    def demag_b(self, b, temperature):
        """return demagnetization Hx, Hy of this element at flux density b
          and temperature"""
        if self.is_magnet():
            # assume polar coordinates of b
            pos = np.arctan2(self.center[1], self.center[0])
            #pos = 0  # cartesian
            mag = self.remanence(temperature)
            magn = np.sqrt(mag[0]**2 + mag[1]**2)
            alfa = np.arctan2(mag[1], mag[0]) - pos
            b1, b2 = b
            bpol = b1 * np.cos(alfa) + b2 * np.sin(alfa)
            reluc = abs(self.reluc[0]) / (4*np.pi*1e-7 * 1000)
            hpol = (bpol - magn)*reluc
            if np.isscalar(hpol):
                if hpol > 0:
                    return 0
            else:
                hpol[hpol > 0] = 0.0
            return -hpol

        return 0

    def permeability(self):
        """return permeability of this element"""
        if self.reluc[0] < 1:
            return 1 / self.reluc[0]
        return 1

    def iron_loss_density(self):
        """return loss_density if element in iron (eg. lamination region)"""
        if self.is_lamination():
            return self.loss_density
        return 0

    def mag_loss_density(self):
        """return loss_density if element in magnet region"""
        if np.any(self.mag):
            return self.loss_density
        return 0

    def wdg_loss_density(self):
        """return loss_density if element in winding region"""
        if self.superelement.subregion:
            if self.superelement.subregion.winding:
                return self.loss_density
        return 0

    def temp(self):
        """return temperature of this element"""
        return sum([v.vpot[1] for v in self.vertices])/len(self.vertices)


class SuperElement(BaseEntity):
    def __init__(self, key, sr_key, elements, nodechains, color,
                 nc_keys, mcvtype, condtype, conduc, length,
                 velsys, velo_1, velo_2, curd_re, curd_im,
                 fillfactor, temp_coef, temperature):
        super(self.__class__, self).__init__(key)
        self.sr_key = sr_key
        self.subregion = None
        self.elements = elements
        for e in elements:
            e.superelement = self
        self.nodechains = nodechains
        self.color = color
        self.nc_keys = nc_keys
        self.mcvtype = mcvtype
        self.condtype = condtype
        self.conduc = conduc
        self.length = length
        self.velsys = velsys
        self.velo = velo_1, velo_2
        self.curd = curd_re, curd_im
        self.fillfactor = fillfactor
        self.temp_coef = temp_coef
        self.temperature = temperature
        self.outside = True

    def area(self):
        """return area of this superelement"""
        return sum([e.area for e in self.elements])

    def get_rect_geom(self):
        """return rectangle parameters of superelement:
        cxy: center coordinates
        w, h: width and height
        alpha: angle of main axis"""
        bxy = np.array([n.xy for b in self.nodechains
                        for n in b.nodes[:-1]])
        # center
        cxy = np.mean(bxy, axis=0)
        # corner points: calculate angles
        b = np.vstack((bxy[-1], bxy, bxy[0]))
        db = np.diff(b, axis=0)
        a = np.arctan2(db[:, 1], db[:, 0])
        da = np.abs(np.diff(a))
        peaks = np.where((da < 6) & (da > 1))[0]
        c = bxy[peaks]
        # width and height (distances between corners)
        dxy = np.linalg.norm(np.vstack(
            (np.diff(bxy, axis=0),
             bxy[-1, :] - bxy[0, :])), axis=1)
        dc = (np.sum(dxy[peaks[0]:peaks[1]]),
              np.sum(dxy[peaks[1]:peaks[2]]),
              np.sum(dxy[peaks[2]:peaks[3]]),
              np.sum(np.hstack(
                  (dxy[peaks[3]:], dxy[:peaks[0]]))))
        w = np.mean(np.sort(dc)[-2:])
        area = self.area()
        h = area/w
        # angle of main axis
        i = np.argmax(dc)
        c = np.vstack((c, c[0]))
        alpha = np.arctan2(c[i+1, 1]-c[i, 1], c[i+1, 0]-c[i, 0])
        if alpha < 0:
            alpha += np.pi
        return {'w': w, 'h': h, 'cxy': cxy,
                'area': area, 'alpha': alpha}


class SubRegion(BaseEntity):
    def __init__(self, key, sr_type, color, name, nturns, curdir, wb_key,
                 superelements, nodechains):
        super(self.__class__, self).__init__(key)
        self.sr_type = sr_type
        self.color = color
        self.name = name
        self.curdir = curdir
        self.num_turns = nturns,
        self.wb_key = wb_key
        self.winding = None
        self.superelements = superelements
        for se in superelements:
            se.subregion = self
        self.nodechains = nodechains

    def elements(self):
        """return elements of this subregion"""
        return [e for s in self.superelements for e in s.elements]

    def outside(self):
        return np.all([se.outside for se in self.superelements])

    def area(self):
        """return area of this subregion"""
        return sum([e.area for e in self.elements()])


class Winding(BaseEntity):
    def __init__(self, key, name, subregions, num_turns, cur_re, cur_im,
                 flux_re, flux_im, volt_re, volt_im):
        super(self.__class__, self).__init__(key)
        self.name = name
        self.subregions = subregions
        for sr in subregions:
            sr.winding = self
        self.num_turns = num_turns
        self.cur = cur_re, cur_im
        self.flux = flux_re, flux_im
        self.volt = volt_re, volt_im

    def elements(self):
        """return elements of this winding"""
        return [e for s in self.subregions for e in s.elements()]


def read(filename):
    """
    Read ISA7 file and return ISA7 object.

    Arguments:
        filename: name of I7/ISA7 file to be read
    """
    import os
    ext = os.path.splitext(filename)[-1]
    if not ext:
        ext = '.I7' if sys.platform == 'win32' else '.ISA7'
        filename += ext
    isa = Isa7(Reader(filename))
    return isa


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    isa = read(filename)
