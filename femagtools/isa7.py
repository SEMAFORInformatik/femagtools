# -*- coding: utf-8 -*-
"""
    femagtools.isa7
    ~~~~~~~~~~~~~~~

    Read FEMAG I7/ISA7 model files
"""
import logging
import struct
import sys
import pdb
import re
import numpy as np
from collections import Counter

logger = logging.getLogger('femagtools.isa7')


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
        self.skip_block(5)
        self.MAGN_TEMPERATURE, self.BR_TEMP_COEF = self.next_block("f")[0:2]
        FC_NUM_CUR_ID, FC_NUM_BETA_ID = self.next_block("i")[0:2]
        if FC_NUM_CUR_ID > 16:
            FC_NUM_CUR_ID = 16

        self.skip_block(3)
        self.skip_block(FC_NUM_CUR_ID * 2)
        self.skip_block(1 + 10 * 5 + 3 + 1 * 5 + 14)

        NUM_FE_EVAL_MOVE_STEP = self.next_block("i")[0]
        if NUM_FE_EVAL_MOVE_STEP < 0:
            NUM_FE_EVAL_MOVE_STEP = 0
            
        self.el_fe_induction_1 = [[[]], [[]], [[]]]
        self.el_fe_induction_2 = [[[]], [[]], [[]]]
        self.eddy_cu_vpot = [[[]], [[]], [[]]]
        self.pos_el_fe_induction = []
        
        if NUM_FE_EVAL_MOVE_STEP > 1:
            self.pos_el_fe_induction = self.next_block("f")
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.el_fe_induction_1[0][0].append(self.next_block("h"))
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
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.eddy_cu_vpot[0][0].append(self.next_block("h"))

        self.skip_block(2)  # start_winkel, end_winkel
        self.skip_block(2 * 5)
        self.skip_block(15)
        self.skip_block(3 * 30 * 30)
        self.skip_block(3)
        self.skip_block(30 * 30)
        self.skip_block(4)
        # stator 3
        self.skip_block(4)
        (yoke_diam, inside_diam,
         slot_height,slot_h1,slot_h2,
         slot_width,slot_r1,slot_r2) = self.next_block("f")[:8]
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
        self.arm_length = arm_length*1e-3  # unit is m
        self.skip_block(2)      
        self.skip_block(30 * 30)
        self.skip_block(30 * 30)
        self.skip_block(1 * 20)
        self.skip_block(8)
        self.beta_loss=self.next_block("h")[:FC_NUM_BETA_ID] # BETA_LOSS_EVAL_STEP
        self.curr_loss=self.next_block("h")[:FC_NUM_CUR_ID] # CURR_LOSS_EVAL_STEP
        FC_NUM_MOVE_LOSSES = self.next_block("i")[0]

        if FC_NUM_MOVE_LOSSES > 1 and NUM_FE_EVAL_MOVE_STEP > 1:
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.el_fe_induction_1[1][0].append(self.next_block("h"))
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
        self.skip_block(14)
        self.skip_block(2 * 3 + 6 * 100 * 3)
        self.skip_block(30)
        self.skip_block(11 * 4)
        self.skip_block()
        self.skip_block(1 * 4)
        # NOM_CURRENT
        # PR_BASIC_LOSS_DATA
        # TOT_MAGNET_AREA
        # MOVE_EXTERN
        # MOVE_ARMATURE
        self.skip_block(5)

        self.pole_pairs, self.poles_sim = self.next_block("i")[:2]
        self.SLOT_WIRE_DIAMETER = self.next_block("f")
        self.SLOT_WIRE_NUMBERS = self.next_block("i")
        self.skip_block(20*(3 + 2 * 20)) # BASE_FREQUENCY ..
        self.skip_block(2) # R_TORQUE .. NUM_NOLOAD_EX_CURRENT_STEPS
        (self.R_CURRENT,
         self.R_LOAD_VOLTAGE,
         self.R_NOLOAD_VOLTAGE) = self.next_block("f")
        x = self.next_block("f")
        self.R_COSPHI = x[0]
        self.R_BETA_OPT = x[1:]
        self.skip_block(10) # R_FLUX_LOAD. NUM_NOLOAD_EX_CURRENT_STEPS

        if (FC_NUM_MOVE_LOSSES > 2 and NUM_FE_EVAL_MOVE_STEP > 1
            and FC_NUM_BETA_ID > 1):
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.el_fe_induction_1[2][0].append(self.next_block("h"))
                self.el_fe_induction_2[2][0].append(self.next_block("h"))
            for i in range(NUM_FE_EVAL_MOVE_STEP + 1):
                self.eddy_cu_vpot[2][0].append(self.next_block("h"))

        self.skip_block()
        self.skip_block(2 * 3)  # MAX_LOSS_EVAL_STEPS
        self.Q_SLOTS_NUMBER, self.M_PHASE_NUMBER = self.next_block("i")[:2]
        self.N_LAYERS_SLOT, self.N_WIRES_PER_SLOT = self.next_block("i")[:2]
        self.skip_block(1)
        self.skip_block(10 * 100) # num_index_cad
        self.skip_block(1 * 100)
        self.skip_block() # index_cad
        self.skip_block(1 * 4) # heat_tranfer_coeff
        self.skip_block(2 * 2)
        self.skip_block()
        self.skip_block(2 * 4)
        self.skip_block(3)
        self.skip_block(1 * 64) # bnodes_mech
        self.skip_block(6)
            
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
        try:
            unpacked = struct.iter_unpack("=" + fmt_,
                                          self.file[self.pos:self.pos
                                                    + blockSize])
            unpacked = [x for x in unpacked]

        except AttributeError:  # python 2 has no iter_unpack
            chunksize = struct.calcsize("=" + fmt_)
            offset = self.pos
            unpacked = []
            for j in range(blockSize // chunksize):
                unpacked.append(struct.unpack_from("=" + fmt_,
                                                   self.file,
                                                   offset))
                offset += chunksize
            logger.info("%s: %d %d", fmt_, blockSize, len(unpacked))

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

    color = {1: [1.0, 0.0, 0.0],
             2: [0.0, 1.0, 0.0],
             3: [1.0, 1.0, 0.0],
             4: [0.0, 0.5019607843137255, 1.0],
             5: [0.9803921568627451, 0.0, 1.0],
             6: [0.0, 1.0, 0.8235294117647058],
             7: [1.0, 1.0, 1.0],
             8: [0.0, 0.0, 0.0],
             9: [0.0, 0.0, 0.5882352941176471],
             10: [0.6666666666666666, 0.0, 0.0],
             11: [0.6666666666666666, 1.0, 0.0],
             12: [1.0, 0.6274509803921569, 0.0],
             13: [0.0, 0.0, 1.0],
             14: [0.6666666666666666, 0.0, 1.0],
             15: [0.0, 0.8235294117647058, 1.0],
             16: [0.8274509803921568, 0.8274509803921568, 0.8274509803921568]}

    def __init__(self, reader):
        self.points = [Point(x, y)
                       for x, y in zip(reader.POINT_ISA_POINT_REC_PT_CO_X,
                                       reader.POINT_ISA_POINT_REC_PT_CO_Y)]

        self.lines = [Line(self.points[abs(pk1) - 1], self.points[abs(pk2) - 1])
                      for pk1, pk2 in zip(reader.LINE_ISA_LINE_REC_LN_PNT_1,
                                          reader.LINE_ISA_LINE_REC_LN_PNT_2)]
        logger.info("Nodes")
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

        logger.info("Nodechains")
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
        logger.info("Elements")
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
            self.elements.append(
                Element(e + 1,
                        reader.ELEM_ISA_ELEM_REC_EL_TYP[e],
                        reader.ELEM_ISA_ELEM_REC_EL_SE_KEY[e] - 1,
                        vertices,
                        (reader.ELEM_ISA_ELEM_REC_EL_RELUC[e],
                         reader.ELEM_ISA_ELEM_REC_EL_RELUC_2[e]),
                        (reader.ELEM_ISA_ELEM_REC_EL_MAG_1[e],
                         reader.ELEM_ISA_ELEM_REC_EL_MAG_2[e]),
                        loss_dens, # in W/m³
                        reader.BR_TEMP_COEF/100)   # in 1/K
            )
        logger.info("SuperElements")
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

            elements = []
            for elk in el_keys:
                elements.append(self.elements[elk - 1])

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
                             reader.SUPEL_ISA_SUPEL_REC_SE_CURD_IM[se]))

        logger.info("Subregions")
        self.subregions = []
        for sr in range(len(reader.SR_ISA_SR_SE_PNTR)):
            se_keys = []
            se_ptr = reader.SR_ISA_SR_SE_PNTR[sr]

            while se_ptr > 0:
                se_keys.append(reader.SR_SE_ISA_SE_KEY[se_ptr - 1])
                se_ptr = reader.SR_SE_ISA_NXT_SE_PNTR[se_ptr - 1]

            superelements = []
            for sek in se_keys:
                superelements.append(self.superelements[sek - 1])

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
                          reader.SR_ISA_SR_REC_SR_NAME[sr],
                          reader.SR_ISA_SR_REC_SR_NTURNS[sr],
                          reader.SR_ISA_SR_REC_SR_CUR_DIR[sr],
                          reader.SR_ISA_SR_REC_SR_WB_KEY[sr] - 1,
                          superelements,
                          nodechains))

        logger.info("Windings")
        self.windings = []
        try:
            for wd in range(len(reader.WB_ISA_WB_SR_PNTR)):
                sr_keys = []
                sr_ptr = reader.WB_ISA_WB_SR_PNTR[wd]

                while sr_ptr > 0:
                    sr_keys.append(reader.WB_SR_ISA_SR_KEY[sr_ptr - 1])
                    sr_ptr = reader.WB_SR_ISA_NXT_SR_PNTR[sr_ptr - 1]

                subregions = []
                for srk in sr_keys:
                    subregions.append(self.subregions[srk - 1])

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

        for a in ('FC_RADIUS', 'pole_pairs', 'poles_sim',
                      'MAGN_TEMPERATURE', 'BR_TEMP_COEF'):
            v = getattr(reader, a, '')
            if v:
                setattr(self, a, v)
        if getattr(reader, 'pole_pairs', 0):
            self.num_poles = 2*self.pole_pairs
        if getattr(reader, 'slots', 0):
            self.num_slots = reader.slots
        try:
            self.arm_length = reader.arm_length*1e-3  # in m
        except:
            pass
        self.pos_el_fe_induction = np.asarray(reader.pos_el_fe_induction)
        try:
            self.beta_loss = np.asarray(reader.beta_loss)
            self.curr_loss = np.array([c/np.sqrt(2) for c in reader.curr_loss]).tolist()
        except AttributeError:
                pass
        if len(np.asarray(reader.el_fe_induction_1).shape) > 2:
            self.el_fe_induction_1 = np.asarray(reader.el_fe_induction_1).T/1000
            self.el_fe_induction_2 = np.asarray(reader.el_fe_induction_2).T/1000
            self.eddy_cu_vpot = np.asarray(reader.eddy_cu_vpot).T/1000
        else:
            self.el_fe_induction_1 = np.asarray(
                [e for e in reader.el_fe_induction_1 if e[0]]).T/1000
            self.el_fe_induction_2 = np.asarray(
                [e for e in reader.el_fe_induction_2 if e[0]]).T/1000
            self.eddy_cu_vpot = np.asarray(
                [e for e in reader.eddy_cu_vpot if e[0]]).T/1000
        logger.info('El Fe Induction %s', np.asarray(reader.el_fe_induction_1).shape)

    def get_subregion(self, name):
        """return subregion by name"""
        for s in self.subregions:
            if s.name == name:
                return s
        raise ValueError('no such subregion "{}" in this model'.format(name))

    def wdg_elements(self):
        """return elements in winding region"""
        return [el for el in self.elements
                if self.superelement.condtype != 0]

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

    def flux_density(self, x, y, icur, ibeta, cosys='cartes'):
        """return pos and flux density (bx, by) or (br, bt)
        at pos x, y for current and beta"""
        el = self.get_element(x, y)
        ekey = el.key-1
        b1 = np.array(self.el_fe_induction_1[ekey, :, icur, ibeta])
        b2 = np.array(self.el_fe_induction_2[ekey, :, icur, ibeta])
        if cosys == 'polar':
            a = np.arctan2(el.center[1], el.center[0])
            br, bphi = np.array(((np.cos(a), -np.sin(a)),
                                     (np.sin(a), np.cos(a)))).dot(
                                         ((b1),(b2)))
            return dict(
                pos = self.pos_el_fe_induction,
                br = br,
                bt = bphi)
        return dict(
            pos = self.pos_el_fe_induction,
            bx = b1,
            by= b2)

    def flux_dens(self, x, y, icur, ibeta, cosys='cartes'):
        return self.flux_density(x, y, icur, ibeta, cosys)
    
    def demagnetization(self, x, y, icur, ibeta, cosys='cartes'):
        el = self.get_element(x, y)
        flxdens = self.flux_density(x,y, icur, ibeta, cosys)
        if cosys == 'polar':
            return (flxdens['pos'], el.demag_b(flxdens['pos'],
                                                   (flxdens['br'], flxdens['bt']),
                                self.MAGN_TEMPERATURE))
        return (flxdens['pos'], el.demag_b(np.arctan2(y, x),
                                               (flxdens['bx'], flxdens['by']),
                                    self.MAGN_TEMPERATURE))
        
    def demag_situation(self, icur, ibeta, hlim):
        """return h max, h avg, area, pos for demag situation for
        each magnet
        Arguments:
          icur: cur amplitude index
          ibeta: beta angle index 
          hlim: limit of demagnetization (kA/m)
        """
        results=[]
        for se in self.magnet_super_elements():
            elements = np.array(se.elements)
            demag = np.array([self.demagnetization(*el.center, icur, ibeta)[1]
                        for el in elements])
            ind = np.unravel_index(np.argmax(demag, axis=None), demag.shape)
            dmax = demag[:, ind[1]]
            area_tot = np.sum([e.area for e in elements])
            area_demag = np.sum([e.area for e in elements[-dmax<hlim]])
            results.append(dict(
                h_max=demag[ind],
                h_avg=np.average(dmax),
                area_tot=area_tot,
                area_demag=area_demag,
                pos=self.pos_el_fe_induction[ind[1]]))
        return results

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
                 se_key, vertices, reluc, mag, loss_density, br_temp_coef=0):
        super(self.__class__, self).__init__(key)
        self.el_type = el_type
        self.se_key = se_key
        self.vertices = vertices
        self.reluc = reluc
        self.mag = mag
        self.br_temp_coef = br_temp_coef
        self.loss_density = loss_density
        if el_type == 1:    # Linear triangle
            self.area = ((vertices[2].x - vertices[1].x) *
                         (vertices[0].y - vertices[1].y) -
                         (vertices[2].y - vertices[1].y) *
                         (vertices[0].x - vertices[1].x))/2
        elif el_type == 2:  # Linear rectangle
            self.area = ((vertices[2].x - vertices[1].x) *
                         (vertices[0].y - vertices[1].y) -
                         (vertices[2].y - vertices[1].y) *
                         (vertices[0].x - vertices[1].x) +
                         (vertices[3].x - vertices[2].x) *
                         (vertices[0].y - vertices[2].y) -
                         (vertices[3].y - vertices[2].y) *
                         (vertices[0].x - vertices[2].x))/2
        elif el_type == 3:  # Square triangle
            self.area = ((vertices[4].x - vertices[2].x) *
                         (vertices[0].y - vertices[2].y) -
                         (vertices[4].y - vertices[1].y) *
                         (vertices[0].x - vertices[2].x))/2
        elif el_type == 4:   # Square rectangle
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
        if self.el_type == 1:
            y31 = ev[2].y - ev[0].y
            y21 = ev[1].y - ev[0].y
            x13 = ev[0].x - ev[2].x
            x21 = ev[1].x - ev[0].x
            a21 = ev[1].vpot[0] - ev[0].vpot[0]
            a31 = ev[2].vpot[0] - ev[0].vpot[0]
            delta = self.superelement.length * (y31 * x21 + y21 * x13)

            b1, b2 = ((x13 * a21 + x21 * a31) / delta,
                    (-y31 * a21 + y21 * a31) / delta)

        elif self.el_type == 2:
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
                                 (-np.sin(a), np.cos(a)))).dot(((b1),(b2)))            
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
    
    def demagnetization(self, temperature=20):
        """return demagnetization of this element"""
        return self.demag_b(0.0, self.flux_density(), temperature)

    def demag_b(self, pos, b, temperature):
        """return demagnetization of this element at flux density b"""
        if self.is_magnet():
            br_temp_corr = 1. +  self.br_temp_coef*(temperature - 20.)
            magn = np.sqrt(self.mag[0]**2 + self.mag[1]**2)*br_temp_corr
            alfa = np.arctan2(self.mag[1], self.mag[0]) - pos
            b1, b2 = b
            bpol = b1 * np.cos(alfa) + b2 * np.sin(alfa)
            reluc = abs(self.reluc[0]) / (4*np.pi*1e-7 * 1000)
            hpol = (bpol - magn)*reluc
            if np.isscalar(hpol):
                if hpol>0:
                    return 0
            else:
                hpol[hpol>0] = 0.0
            return -hpol
            
        return 0

    def permeability(self):
        """return permeability of this element"""
        if self.reluc[0] < 1:
            return 1 / self.reluc[0]
        return 1

    def iron_loss_density(self):
        """return loss_density if element in iron (eg. lamination region)"""
        if self.reluc != (1.0, 1.0) and self.mag == (0.0, 0.0):
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
        

class SuperElement(BaseEntity):
    def __init__(self, key, sr_key, elements, nodechains, color,
                 nc_keys, mcvtype, condtype, conduc, length,
                 velsys, velo_1, velo_2, curd_re, curd_im):
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
        return [e for s in self.subregions for e in s.elements]


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
