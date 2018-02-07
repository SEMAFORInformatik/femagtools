""" femagtools.isa7 """
import logging
import struct
import sys
import numpy as np
from collections import Counter

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
    
    def __init__(self, filename):
        self.read(filename)
        
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
            
            point1 = self.points[abs(pk1)-1]
            point2 = self.points[abs(pk2)-1]
            
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

            node1 = self.nodes[abs(nd1)-1]
            nodem = self.nodes[ndm-1]
            node2 = self.nodes[abs(nd2)-1]

            if nd1 < 0 or nd2 < 0:
                nodes = node1, nodem, node2
            elif ndm > 0:
                nodes = node1, nodem, node2
            else:
                nodes = node1, None, node2
                
            self.nodechains.append(
                NodeChain(self.NDCHN_ISA_NC_VALID, nc+1, nodes))
            
        self.elements = []
        for e in range(self.NUM_ELE):
            ndkeys = []
            ndk = self.ELEM_ISA_EL_NOD_PNTR[e]

            while ndk > 0:
                ndkeys.append(self.ELE_NOD_ISA_ND_KEY[ndk-1])
                ndk = self.ELE_NOD_ISA_NXT_ND_PNTR[ndk-1]

            vertices = [self.nodes[k-1] for k in ndkeys + [ndkeys[0]]]

            self.elements.append(
                Element(self.ELEM_ISA_EL_VALID[e],
                        self.ELEM_ISA_ELEM_REC_EL_TYP[e],
                        self.ELEM_ISA_ELEM_REC_EL_SE_KEY[e]-1,
                        vertices))

        self.superelements = []
        for se in range(self.NUM_SPEL):
            nc_keys = []
            nc_ptr = self.SUPEL_ISA_SE_NDCHN_PNTR[se]

            while nc_ptr > 0:
                nc_keys.append(self.SE_NDCHN_ISA_NC_KEY[nc_ptr-1])
                nc_ptr = self.SE_NDCHN_ISA_NXT_NC_PNTR[nc_ptr-1]
                
            nodechains = []
            for nck in nc_keys:
                if nck > 0:
                    nodechains.append(self.nodechains[abs(nck)-1])
                else:
                    nodechains.append(self.nodechains[abs(nck)-1].reverse())
                    
            el_keys = []
            el_ptr = self.SUPEL_ISA_SE_EL_PNTR[se]
            
            while el_ptr > 0:
                el_keys.append(self.SE_EL_ISA_EL_KEY[el_ptr-1])
                el_ptr = self.SE_EL_ISA_NXT_EL_PNTR[el_ptr-1]
                
            elements = []
            for elk in el_keys:
                elements.append(self.elements[elk-1])
            
            self.superelements.append(
                SuperElement(self.SUPEL_ISA_SE_VALID[se],
                             self.SUPEL_ISA_SUPEL_REC_SE_SR_KEY[se]-1,
                             elements,
                             nodechains,
                             self.SUPEL_ISA_SUPEL_REC_SE_COL[se],
                             nc_keys))
            
        self.subregions = []
        for sr in range(self.NUM_SR):
            se_keys = []
            se_ptr = self.SR_ISA_SR_SE_PNTR[sr]
            
            while se_ptr > 0:
                se_keys.append(self.SR_SE_ISA_SE_KEY[se_ptr-1])
                se_ptr = self.SR_SE_ISA_NXT_SE_PNTR[se_ptr-1]
                
            superelements = []
            for sek in se_keys:
                superelements.append(self.superelements[sek-1])
                
            nodechains = []
            nc_keys = []
            for se in superelements:
                nc_keys.extend([abs(nc.key) for nc in se.nodechains])
            nc_keys = [nck for nck, count in Counter(nc_keys).items() if count < 2]
            for se in superelements:
                nodechains.extend([nc for nc in se.nodechains if abs(nc.key) in nc_keys])
                
            self.subregions.append(
                SubRegion(self.SR_ISA_SR_VALID[sr],
                          self.SR_ISA_SR_REC_SR_TYP[sr],
                          self.SR_ISA_SR_REC_SR_COL[sr],
                          self.SR_ISA_SR_REC_SR_NAME[sr],
                          self.SR_ISA_SR_REC_SR_CUR_DIR[sr],
                          self.SR_ISA_SR_REC_SR_WB_KEY[sr]-1,
                          superelements,
                          nodechains))
            
            
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
         self.SR_ISA_SR_REC_SR_VOLT_IM) = reader.next_block("?hhhihhhhfffff")
        
        (self.SR_SE_ISA_SE_KEY,
         self.SR_SE_ISA_NXT_SE_PNTR) = reader.next_block("hh")


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
        
        
class Node(object):
    def __init__(self, valid, bndcnd, pernod, x, y, vpot_re, vpot_im):
        self.valid = valid
        self.bndcnd = bndcnd,
        self.pernod = pernod,
        self.x = x
        self.y = y
        self.xy = x, y
        self.vpot = vpot_re, vpot_im
        

class NodeChain(object):
    def __init__(self, valid, key, nodes):
        self.valid = valid
        self.key = key
        self.node1 = nodes[0]
        self.nodemid = nodes[1]   
        self.node2 = nodes[2]
        self.nodes = (nodes[0], nodes[2]) if nodes[1] is None else (nodes[0], nodes[1], nodes[2])
    
    def reverse(self):
        return NodeChain(self.valid, self.key * (-1),
                         [self.node2, self.nodemid, self.node1])
    
    
class Element(object):
    def __init__(self, valid, el_type, se_key, vertices):
        self.valid = valid
        self.el_type = el_type
        self.se_key = se_key
        self.vertices = vertices
        
        
class SuperElement(object):
    def __init__(self, valid, sr_key, elements, nodechains,
                 color, nc_keys):
        self.valid = valid
        self.sr_key = sr_key
        self.elements = elements
        self.nodechains = nodechains
        self.color = color
        self.nc_keys = nc_keys
        
        
class SubRegion(object):
    def __init__(self, valid, sr_type, color, name, curdir,
                 wb_key, superelements, nodechains):
        self.valid = valid
        self.sr_type = sr_type
        self.color = color
        self.name = name
        self.curdir = curdir
        self.wb_key = wb_key
        self.superelements = superelements
        self.nodechains = nodechains

        
def read(filename):
    """Read ISA7 file and return ISA7 object."""
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
    
