# -*- coding: utf-8 -*-
"""
    femagtools.isa7
    ~~~~~~~~~~~~~~~

    Read FEMAG I7/ISA7 model files
"""
import logging
import struct
import sys
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
        with open(filename, mode="rb") as self.file:
            self.file = self.file.read()
        self.pos = 0

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
        for s in re.findall("[0-9]*.|[0-9]*\?", fmt):
            if len(s) > 1 and s[-1] != "s":
                fmt_ += int(s[:-1]) * s[-1]
            else:
                fmt_ += s
        values = []
        for i, dtype in enumerate(re.findall("\?|[0-9]*s?", fmt_)[:-1]):
            if dtype == "?":
                values.append([bool(u[i]) for u in unpacked])
            elif "s" in dtype:
                values.append([u[i].decode() for u in unpacked])
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

    def __init__(self, filename):
        logger.info("read %s", filename)
        self.__read(filename)

        self.points = []
        for p in range(self.NUM_PNT):
            self.points.append(
                Point(self.POINT_ISA_PT_VALID,
                      self.POINT_ISA_POINT_REC_PT_CO_X[p],
                      self.POINT_ISA_POINT_REC_PT_CO_Y[p]))

        self.lines = []
        for ln in range(self.NUM_LIN):
            pk1 = self.LINE_ISA_LINE_REC_LN_PNT_1[ln]
            pk2 = self.LINE_ISA_LINE_REC_LN_PNT_2[ln]

            point1 = self.points[abs(pk1) - 1]
            point2 = self.points[abs(pk2) - 1]

            if pk1 > 0 and pk2 > 0:
                self.lines.append(
                    Line(self.LINE_ISA_LN_VALID, point1, point2))
            else:
                self.lines.append(
                    Line(self.LINE_ISA_LN_VALID, point1, point2))

        self.nodes = []
        for n in range(self.NUM_NOD):
            self.nodes.append(
                Node(self.NODE_ISA_ND_VALID[n],
                     n + 1,
                     self.NODE_ISA_NODE_REC_ND_BND_CND[n],
                     self.NODE_ISA_NODE_REC_ND_PER_NOD[n],
                     self.NODE_ISA_NODE_REC_ND_CO_1[n],
                     self.NODE_ISA_NODE_REC_ND_CO_2[n],
                     self.NODE_ISA_NODE_REC_ND_VP_RE[n],
                     self.NODE_ISA_NODE_REC_ND_VP_IM[n]))

        self.nodechains = []
        for nc in range(self.NUM_NDCH):
            nd1 = self.NDCHN_ISA_NDCHN_REC_NC_NOD_1[nc]
            nd2 = self.NDCHN_ISA_NDCHN_REC_NC_NOD_2[nc]
            ndm = self.NDCHN_ISA_NDCHN_REC_NC_NOD_MID[nc]

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
                NodeChain(self.NDCHN_ISA_NC_VALID, nc + 1, nodes))

        self.elements = []
        for e in range(self.NUM_ELE):
            ndkeys = []
            ndk = self.ELEM_ISA_EL_NOD_PNTR[e]

            while ndk > 0:
                ndkeys.append(self.ELE_NOD_ISA_ND_KEY[ndk - 1])
                ndk = self.ELE_NOD_ISA_NXT_ND_PNTR[ndk - 1]

            vertices = [self.nodes[k - 1] for k in ndkeys]

            self.elements.append(
                Element(self.ELEM_ISA_EL_VALID[e],
                        e + 1,
                        self.ELEM_ISA_ELEM_REC_EL_TYP[e],
                        self.ELEM_ISA_ELEM_REC_EL_SE_KEY[e] - 1,
                        vertices,
                        (self.ELEM_ISA_ELEM_REC_EL_RELUC[e],
                         self.ELEM_ISA_ELEM_REC_EL_RELUC_2[e]),
                        (self.ELEM_ISA_ELEM_REC_EL_MAG_1[e],
                         self.ELEM_ISA_ELEM_REC_EL_MAG_2[e]),
                        self.ELEM_ISA_ELEM_REC_LOSS_DENS[e] * 1e-6)
            )

        self.superelements = []
        for se in range(self.NUM_SPEL):
            nc_keys = []
            nc_ptr = self.SUPEL_ISA_SE_NDCHN_PNTR[se]

            while nc_ptr > 0:
                nc_keys.append(self.SE_NDCHN_ISA_NC_KEY[nc_ptr - 1])
                nc_ptr = self.SE_NDCHN_ISA_NXT_NC_PNTR[nc_ptr - 1]

            nodechains = []
            for nck in nc_keys:
                if nck > 0:
                    nodechains.append(self.nodechains[abs(nck) - 1])
                else:
                    nodechains.append(self.nodechains[abs(nck) - 1].reverse())

            el_keys = []
            el_ptr = self.SUPEL_ISA_SE_EL_PNTR[se]

            while el_ptr > 0:
                el_keys.append(self.SE_EL_ISA_EL_KEY[el_ptr - 1])
                el_ptr = self.SE_EL_ISA_NXT_EL_PNTR[el_ptr - 1]

            elements = []
            for elk in el_keys:
                elements.append(self.elements[elk - 1])

            self.superelements.append(
                SuperElement(self.SUPEL_ISA_SE_VALID[se],
                             se + 1,
                             self.SUPEL_ISA_SUPEL_REC_SE_SR_KEY[se] - 1,
                             elements,
                             nodechains,
                             self.SUPEL_ISA_SUPEL_REC_SE_COL[se],
                             nc_keys,
                             self.SUPEL_ISA_SUPEL_REC_SE_MCV_TYP[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_COND_TYP[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_CONDUC[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_LENGHT[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_VEL_SYS[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_VELO_1[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_VELO_2[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_CURD_RE[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_CURD_IM[se]))

        self.subregions = []
        for sr in range(self.NUM_SR):
            se_keys = []
            se_ptr = self.SR_ISA_SR_SE_PNTR[sr]

            while se_ptr > 0:
                se_keys.append(self.SR_SE_ISA_SE_KEY[se_ptr - 1])
                se_ptr = self.SR_SE_ISA_NXT_SE_PNTR[se_ptr - 1]

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
                SubRegion(self.SR_ISA_SR_VALID[sr],
                          sr + 1,
                          self.SR_ISA_SR_REC_SR_TYP[sr],
                          self.SR_ISA_SR_REC_SR_COL[sr],
                          self.SR_ISA_SR_REC_SR_NAME[sr],
                          self.SR_ISA_SR_REC_SR_CUR_DIR[sr],
                          self.SR_ISA_SR_REC_SR_WB_KEY[sr] - 1,
                          superelements,
                          nodechains))

        self.windings = []
        for wd in range(self.NUM_WB):
            sr_keys = []
            sr_ptr = self.WB_ISA_WB_SR_PNTR[wd]

            while sr_ptr > 0:
                sr_keys.append(self.WB_SR_ISA_SR_KEY[sr_ptr - 1])
                sr_ptr = self.WB_SR_ISA_NXT_SR_PNTR[sr_ptr - 1]

            subregions = []
            for srk in sr_keys:
                subregions.append(self.subregions[srk - 1])

            self.windings.append(
                Winding(self.WB_ISA_WB_VALID[wd],
                        wd + 1,
                        self.WB_ISA_WB_REC_WB_NAME[wd],
                        subregions,
                        self.WB_ISA_WB_REC_WB_TURN[wd],
                        self.WB_ISA_WB_REC_WB_GCUR_RE[wd],
                        self.WB_ISA_WB_REC_WB_GCUR_IM[wd],
                        self.WB_ISA_WB_REC_WB_IMPDZ_RE[wd],
                        self.WB_ISA_WB_REC_WB_IMPDZ_IM[wd],
                        self.WB_ISA_WB_REC_WB_VOLT_RE[wd],
                        self.WB_ISA_WB_REC_WB_VOLT_IM[wd]))
        logger.info("Total nodes %d elements %d superelements %d subregions %d",
                    len(self.nodes), len(self.elements),
                    len(self.superelements),
                    len(self.subregions))

    def __read(self, filename):
        reader = Reader(filename)

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
         ) = reader.next_block("i")

        (self.POINT_ISA_PT_VALID,
         self.POINT_ISA_POINT_REC_PT_CO_X,
         self.POINT_ISA_POINT_REC_PT_CO_Y) = reader.next_block("?ff")

        (self.LINE_ISA_LN_VALID,
         self.LINE_ISA_LINE_REC_LN_PNT_1,
         self.LINE_ISA_LINE_REC_LN_PNT_2) = reader.next_block("?hh")

        (self.NODE_ISA_ND_VALID,
         self.NODE_ISA_NOD_EL_PNTR,
         self.NODE_ISA_ND_CO_RAD,
         self.NODE_ISA_ND_CO_PHI,
         self.NODE_ISA_NODE_REC_ND_BND_CND,
         self.NODE_ISA_NODE_REC_ND_PER_NOD,
         self.NODE_ISA_NODE_REC_ND_SV_PNTR,
         self.NODE_ISA_NODE_REC_ND_CO_1,
         self.NODE_ISA_NODE_REC_ND_CO_2,
         self.NODE_ISA_NODE_REC_ND_VP_RE,
         self.NODE_ISA_NODE_REC_ND_VP_IM) = reader.next_block("?iffhiiffff")

        (self.NOD_ELE_ISA_EL_KEY,
         self.NOD_ELE_ISA_NXT_EL_PNTR) = reader.next_block("ii")

        (self.NDCHN_ISA_NC_VALID,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_1,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_2,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_MID) = reader.next_block("?iii")

        (self.ELEM_ISA_EL_VALID,
         self.ELEM_ISA_EL_NOD_PNTR,
         self.ELEM_ISA_ELEM_REC_EL_TYP,
         self.ELEM_ISA_ELEM_REC_EL_SE_KEY,
         self.ELEM_ISA_ELEM_REC_EL_RELUC,
         self.ELEM_ISA_ELEM_REC_EL_RELUC_2,
         self.ELEM_ISA_ELEM_REC_EL_MAG_1,
         self.ELEM_ISA_ELEM_REC_EL_MAG_2) = reader.next_block("?ihhffff")

        (self.ELE_NOD_ISA_ND_KEY,
         self.ELE_NOD_ISA_NXT_ND_PNTR) = reader.next_block("ii")

        (self.SUPEL_ISA_SE_VALID,
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
         ) = reader.next_block("?iihhhhhffffff")

        (self.SE_NDCHN_ISA_NC_KEY,
         self.SE_NDCHN_ISA_NXT_NC_PNTR) = reader.next_block("ii")

        (self.SE_EL_ISA_EL_KEY,
         self.SE_EL_ISA_NXT_EL_PNTR) = reader.next_block("ii")

        (self.SR_ISA_SR_VALID,
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
         self.SR_ISA_SR_REC_SR_VOLT_IM) = reader.next_block("?hhh4shhhhfffff")

        (self.SR_SE_ISA_SE_KEY,
         self.SR_SE_ISA_NXT_SE_PNTR) = reader.next_block("hh")

        (self.WB_ISA_WB_VALID,
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
         self.WB_ISA_WB_REC_WB_IMPDZ_IM) = reader.next_block("?hh4shhhfffffff")

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
         self.WB_SR_ISA_NXT_SR_PNTR) = reader.next_block("hh")

        reader.skip_block(21)

        ANZAHL_TG = reader.next_block("i")[1]

        reader.skip_block(7)
        reader.skip_block(ANZAHL_TG + 1)
        reader.skip_block(1)

        self.FC_RADIUS = reader.next_block("f")[0]

        reader.skip_block(9)
        FC_NUM_CUR_ID, FC_NUM_BETA_ID = reader.next_block("i")[0:2]
        if FC_NUM_CUR_ID > 16:
            FC_NUM_CUR_ID = 16

        reader.skip_block(3)
        reader.skip_block(FC_NUM_CUR_ID * 2)
        reader.skip_block(1 + 10 * 5 + 3 + 1 * 5 + 14)

        NUM_FE_EVAL_MOVE_STEP = reader.next_block("i")[0]
        if NUM_FE_EVAL_MOVE_STEP < 0:
            NUM_FE_EVAL_MOVE_STEP = 0

        if NUM_FE_EVAL_MOVE_STEP > 1:
            reader.skip_block()
            reader.skip_block((NUM_FE_EVAL_MOVE_STEP + 1) * 2)

        FC_NUM_MOVE_CALC_LOAD_PMS, FC_NUM_FLX = reader.next_block("i")[0:2]

        if FC_NUM_MOVE_CALC_LOAD_PMS > 1:
            reader.skip_block(4)
            reader.skip_block(3 * FC_NUM_FLX)
            reader.skip_block()

        FC_NUM_MOVE_NOLOAD_PMS = reader.next_block("i")[0]

        if FC_NUM_MOVE_NOLOAD_PMS > 1:
            reader.skip_block(4)
            reader.skip_block(2 * FC_NUM_FLX)
            reader.skip_block()

        if NUM_FE_EVAL_MOVE_STEP > 1:
            reader.skip_block(NUM_FE_EVAL_MOVE_STEP + 1)

        reader.skip_block(2)
        reader.skip_block(2 * 5)
        reader.skip_block(15)
        reader.skip_block(3 * 30 * 30)
        reader.skip_block(3)
        reader.skip_block(30 * 30)
        reader.skip_block(21)
        reader.skip_block(30 * 30)
        reader.skip_block(30 * 30)
        reader.skip_block(1 * 20)
        reader.skip_block(10)

        FC_NUM_MOVE_LOSSES = reader.next_block("i")[0]

        if FC_NUM_MOVE_LOSSES > 1 and NUM_FE_EVAL_MOVE_STEP > 1:
            reader.skip_block(2 * (NUM_FE_EVAL_MOVE_STEP + 1))
            reader.skip_block(NUM_FE_EVAL_MOVE_STEP + 1)

        reader.skip_block(74)

        ANZ_FORCE_AREAS = reader.next_block("i")[0]

        if ANZ_FORCE_AREAS > 3:
            ANZ_FORCE_AREAS = 3

        reader.skip_block()
        reader.skip_block(2 * ANZ_FORCE_AREAS)
        reader.skip_block(14)
        reader.skip_block(2 * 3 + 6 * 100 * 3)
        reader.skip_block(30)
        reader.skip_block(11 * 4)
        reader.skip_block()
        reader.skip_block(1 * 4)
        reader.skip_block(8)
        reader.skip_block(3 * 20 + 2 * 20 * 20)
        reader.skip_block(14)

        if (FC_NUM_MOVE_LOSSES > 2 and NUM_FE_EVAL_MOVE_STEP > 1
           and FC_NUM_BETA_ID > 2):
            reader.skip_block(2 * NUM_FE_EVAL_MOVE_STEP + 1)
            reader.skip_block(1 * NUM_FE_EVAL_MOVE_STEP + 1)
            reader.skip_block()
            reader.skip_block(1 * NUM_FE_EVAL_MOVE_STEP + 1)

        reader.skip_block()
        reader.skip_block(2 * 3)
        reader.skip_block(3)
        reader.skip_block(10 * 100)
        reader.skip_block(1 * 100)
        reader.skip_block()
        reader.skip_block(1 * 4)
        reader.skip_block(2 * 2)
        reader.skip_block()
        reader.skip_block(2 * 4)
        reader.skip_block(3)
        reader.skip_block(1 * 64)
        reader.skip_block(6)

        self.ELEM_ISA_ELEM_REC_LOSS_DENS = reader.next_block("f")


class Point(object):
    def __init__(self, valid, x, y):
        self.valid = valid
        self.x = x
        self.y = y
        self.xy = x, y


class Line(object):
    def __init__(self, valid, p1, p2):
        self.valid = valid
        self.p1 = p1
        self.p2 = p2


class BaseEntity(object):
    def __init__(self, valid, key):
        self.valid = valid
        self.key = key


class Node(BaseEntity):
    def __init__(self, valid, key, bndcnd, pernod, x, y, vpot_re, vpot_im):
        super(self.__class__, self).__init__(valid, key)
        self.bndcnd = bndcnd
        self.pernod = pernod
        self.x = x
        self.y = y
        self.xy = x, y
        self.vpot = vpot_re, vpot_im

    def on_boundary(self):
        return self.bndcnd != 0 or self.pernod != 0


class NodeChain(BaseEntity):
    def __init__(self, valid, key, nodes):
        super(self.__class__, self).__init__(valid, key)
        self.node1 = nodes[0]
        self.nodemid = nodes[1]
        self.node2 = nodes[2]
        if nodes[1] is None:
            self.nodes = (nodes[0], nodes[2])
        else:
            self.nodes = (nodes[0], nodes[1], nodes[2])

    def reverse(self):
        return NodeChain(self.valid, self.key * (-1),
                         [self.node2, self.nodemid, self.node1])


class Element(BaseEntity):
    def __init__(self, valid, key, el_type,
                 se_key, vertices, reluc, mag, loss_density):
        super(self.__class__, self).__init__(valid, key)
        self.el_type = el_type
        self.se_key = se_key
        self.vertices = vertices
        self.reluc = reluc
        self.mag = mag
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

    def induction(self):
        """return induction components of this element"""
        ev = self.vertices
        if self.el_type == 1:
            y31 = ev[2].y - ev[0].y
            y21 = ev[1].y - ev[0].y
            x13 = ev[0].x - ev[2].x
            x21 = ev[1].x - ev[0].x
            a21 = ev[1].vpot[0] - ev[0].vpot[0]
            a31 = ev[2].vpot[0] - ev[0].vpot[0]
            delta = self.se_length * (y31 * x21 + y21 * x13)

            return ((x13 * a21 + x21 * a31) / delta,
                    (-y31 * a21 + y21 * a31) / delta)

        elif self.el_type == 2:
            y31 = ev[2].y - ev[0].y
            y21 = ev[1].y - ev[0].y
            x13 = ev[0].x - ev[2].x
            x21 = ev[1].x - ev[0].x
            a21 = ev[1].vpot[0] - ev[0].vpot[0]
            a31 = ev[2].vpot[0] - ev[0].vpot[0]
            delta = self.se_length * (y31 * x21 + y21 * x13)
            b1_a = (x13 * a21 + x21 * a31) / delta
            b2_a = (y21 * a31 - y31 * a21) / delta

            y31 = ev[0].y - ev[2].y
            y21 = ev[3].y - ev[2].y
            x13 = ev[2].x - ev[0].x
            x21 = ev[3].x - ev[2].x
            a24 = ev[3].vpot[0] - ev[2].vpot[0]
            a34 = ev[0].vpot[0] - ev[2].vpot[0]
            delta = self.se_length * (y31 * x21 + y21 * x13)
            b1_b = (x13 * a24 + x21 * a34) / delta
            b2_b = (y21 * a34 - y31 * a24) / delta

            return ((b1_a + b1_b) / 2,
                    (b2_a + b2_b) / 2)

        return (0, 0)

    def demagnetization(self):
        """return demagnetization of this element"""
        if abs(self.mag[0]) > 1e-5 or abs(self.mag[1]) > 1e-5:
            magn = np.sqrt(self.mag[0]**2 + self.mag[1]**2)
            alfa = np.arctan2(self.mag[1], self.mag[0])
            b1, b2 = self.induction()
            bpol = b1 * np.cos(alfa) + b2 * np.sin(alfa)
            hpol = bpol - magn
            if hpol < 0:
                reluc = abs(self.reluc[0]) / (4*np.pi*1e-7 * 1000)
                return abs(hpol * reluc)
        return 0

    def permeability(self):
        """return permeability of this element"""
        if self.reluc[0] < 1:
            return 1 / self.reluc[0]
        return 0

    def loss_density(self):
        return self.loss_density


class SuperElement(BaseEntity):
    def __init__(self, valid, key, sr_key, elements, nodechains, color,
                 nc_keys, mcvtype, condtype, conduc, length,
                 velsys, velo_1, velo_2, curd_re, curd_im):
        super(self.__class__, self).__init__(valid, key)
        self.sr_key = sr_key
        self.elements = elements
        for e in elements:
            e.se_length = length
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
    def __init__(self, valid, key, sr_type, color, name, curdir, wb_key,
                 superelements, nodechains):
        super(self.__class__, self).__init__(valid, key)
        self.sr_type = sr_type
        self.color = color
        self.name = name
        self.curdir = curdir
        self.wb_key = wb_key
        self.superelements = superelements
        self.nodechains = nodechains


class Winding(BaseEntity):
    def __init__(self, valid, key, name, subregions, num_turns, cur_re, cur_im,
                 flux_re, flux_im, volt_re, volt_im):
        super(self.__class__, self).__init__(valid, key)
        self.name = name
        self.subregions = subregions
        self.num_turns = num_turns
        self.cur = cur_re, cur_im
        self.flux = flux_re, flux_im
        self.volt = volt_re, volt_im


def read(filename):
    """
    Read ISA7 file and return ISA7 object.

    Arguments:
        filename: name of I7/ISA7 file to be read
    """
    isa = Isa7(filename)
    return isa


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    isa = read(filename)
