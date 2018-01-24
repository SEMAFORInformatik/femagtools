""" femagtools.isa7 """
import logging
import struct
import sys

logger = logging.getLogger('femagtools.isa7')


class Reader(object):
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
                                          self.file[self.pos:self.pos +
                                                    blockSize])
            unpacked = [x for x in unpacked]

        except AttributeError:  # python 2 has no iter_unpack
            chunksize = struct.calcsize("=" + fmt_)
            offset = self.pos
            unpacked = []
            for j in range(blockSize//chunksize):
                unpacked.append(struct.unpack_from("=" + fmt_, self.file, offset))
                offset += chunksize
            logger.info("%s: %d %d", fmt_, blockSize, len(unpacked))
            
        values = [[bool(x[i]) if fmt[i] == "?" else x[i]
                   for x in unpacked] for i in range(len(fmt))]
        self.pos += blockSize + 4
        if len(fmt) == 1:
            return values[0]
        return values


class Isa7(object):
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
    
    def __init__(self):
        pass
    
    def read(self, filename):
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
         self.NUM_WN_SW, self.WN_SW_PTR, self.WN_SW_HIDX) = reader.next_block("i")

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
         self.SUPEL_ISA_SUPEL_REC_SE_CURD_IM) = reader.next_block("?iihhhhhffffff")

        (self.SE_NDCHN_ISA_NC_KEY,
         self.SE_NDCHN_ISA_NXT_NC_PNTR) = reader.next_block("ii")

        (self.SE_EL_ISA_EL_KEY,
         self.SE_EL_ISA_NXT_EL_PNTR) = reader.next_block("ii")

    def lines(self):
        """Return list of lines"""
        ln = []
        for i in range(self.NUM_LIN):
            x1 = self.POINT_ISA_POINT_REC_PT_CO_X[
                abs(self.LINE_ISA_LINE_REC_LN_PNT_1[i])-1]
            x2 = self.POINT_ISA_POINT_REC_PT_CO_X[
                abs(self.LINE_ISA_LINE_REC_LN_PNT_2[i])-1]
            y1 = self.POINT_ISA_POINT_REC_PT_CO_Y[
                abs(self.LINE_ISA_LINE_REC_LN_PNT_1[i])-1]
            y2 = self.POINT_ISA_POINT_REC_PT_CO_Y[
                abs(self.LINE_ISA_LINE_REC_LN_PNT_2[i])-1]
            p1 = (x1, y1)
            p2 = (x2, y2)
            ln.append((p1, p2))
        return ln

    def mesh(self):
        """Return element outlines"""
        segments = []
        for e in range(self.NUM_ELE):
            ndkeys = []
            nd = 0
            i = self.ELEM_ISA_EL_NOD_PNTR[e]
            
            while i > 0 and nd <= 9:
                nd += 1
                ndkeys.append(self.ELE_NOD_ISA_ND_KEY[i-1])
                i = self.ELE_NOD_ISA_NXT_ND_PNTR[i-1]
    
                segments.append([[self.NODE_ISA_NODE_REC_ND_CO_1[j-1],
                                  self.NODE_ISA_NODE_REC_ND_CO_2[j-1]]
                                 for j in ndkeys+[ndkeys[0]]])
        return segments
    
    def se_outline(self, spels=None):
        """Return list of nodechains"""
        outln = []
        if spels is None:
            spels = range(self.NUM_SPEL)
        elif type(spels) == int:
            spels = [spels]

        for se in spels:
            ncp = self.SUPEL_ISA_SE_NDCHN_PNTR[se]
            nckeys = []

            while ncp > 0:
                nckeys.append(self.SE_NDCHN_ISA_NC_KEY[ncp-1])
                ncp = self.SE_NDCHN_ISA_NXT_NC_PNTR[ncp-1]

            for nck in nckeys:
                if nck > 0:
                    nd1 = abs(self.NDCHN_ISA_NDCHN_REC_NC_NOD_1[abs(nck)-1])
                    nd2 = abs(self.NDCHN_ISA_NDCHN_REC_NC_NOD_2[abs(nck)-1])
                elif nck < 0:
                    nd1 = abs(self.NDCHN_ISA_NDCHN_REC_NC_NOD_2[abs(nck)-1])
                    nd2 = abs(self.NDCHN_ISA_NDCHN_REC_NC_NOD_1[abs(nck)-1])

                p1 = (self.NODE_ISA_NODE_REC_ND_CO_1[nd1-1],
                      self.NODE_ISA_NODE_REC_ND_CO_2[nd1-1])
                p2 = (self.NODE_ISA_NODE_REC_ND_CO_1[nd2-1],
                      self.NODE_ISA_NODE_REC_ND_CO_2[nd2-1])

                outln.append((p1, p2))

        return outln


def read(filename):
    """Read ISA7 file and return ISA7 object."""
    isa = Isa7()
    isa.read(filename)
    return isa

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    isa = read(filename)
    
