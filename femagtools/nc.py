# -*- coding: utf-8 -*-
"""
    femagtools.nc
    ~~~~~~~~~~~~~

    Read FEMAG netcdf model files
"""
import logging
import netCDF4
from femagtools import isa7


logger = logging.getLogger('femagtools.nc')


class Reader(object):
    """
    Open and Read NetCDF file

    Arguments:
        filename: name of NetCDF nc file to be read
    """

    def __init__(self, filename):
        ds = netCDF4.Dataset(filename)
        self.POINT_ISA_POINT_REC_PT_CO_X = []
        self.POINT_ISA_POINT_REC_PT_CO_Y = []
        self.LINE_ISA_LINE_REC_LN_PNT_1 = []
        self.LINE_ISA_LINE_REC_LN_PNT_2 = []
        try:
            grp = ds.groups['points']
            (self.POINT_ISA_POINT_REC_PT_CO_X,
             self.POINT_ISA_POINT_REC_PT_CO_Y) = [
                 grp.variables[k][:] for k in ('co_x', 'co_y')]
        except KeyError:
            pass
        try:
            grp = ds.groups['lines']
            (self.LINE_ISA_LINE_REC_LN_PNT_1,
             self.LINE_ISA_LINE_REC_LN_PNT_2) = [
                 grp.variables[k][:] for k in ('pnt_1', 'pnt_2')]
        except KeyError:
            pass
        grp = ds.groups['nodes']
        (self.NODE_ISA_NOD_EL_PNTR,
         self.NODE_ISA_NODE_REC_ND_BND_CND,
         self.NODE_ISA_NODE_REC_ND_PER_NOD,
         self.NODE_ISA_NODE_REC_ND_CO_1,
         self.NODE_ISA_NODE_REC_ND_CO_2,
         self.NODE_ISA_ND_CO_RAD,
         self.NODE_ISA_ND_CO_PHI,
         self.NODE_ISA_NODE_REC_ND_VP_RE,
         self.NODE_ISA_NODE_REC_ND_VP_IM) = [
             grp.variables[k][:-1]
             for k in ('bnd_cnd', 'bnd_cnd', 'per_nod',
                       'co_1', 'co_2', 'co_rad', 'co_phi', 'vp_re', 'vp_im')]
        logger.debug('Nodes: %d', len(self.NODE_ISA_NODE_REC_ND_CO_1))
        
        grp = ds.groups['node_elements']
        (self.NODE_ELE_ISA_NOD_EL_KEY,
         self.NODE_ELE_ISA_NOD_NXT_EL_PNTR) = [
             grp.variables[k][:-1]
             for k in ('el_key', 'nxt_el_pntr')]

        grp = ds.groups['nodechains']
        (self.NDCHN_ISA_NDCHN_REC_NC_NOD_1,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_2,
         self.NDCHN_ISA_NDCHN_REC_NC_NOD_MID) = [
             grp.variables[k][:-1]
             for k in ('nod_1', 'nod_2', 'nod_mid')]

        grp = ds.groups['elements']
        (self.ELEM_ISA_EL_NOD_PNTR,
         self.ELEM_ISA_ELEM_REC_EL_TYP,
         self.ELEM_ISA_ELEM_REC_EL_SE_KEY,
         self.ELEM_ISA_ELEM_REC_EL_RELUC,
         self.ELEM_ISA_ELEM_REC_EL_RELUC_2,
         self.ELEM_ISA_ELEM_REC_EL_MAG_1,
         self.ELEM_ISA_ELEM_REC_EL_MAG_2,
         self.ELEM_ISA_ELEM_REC_LOSS_DENS) = [
             grp.variables[k][:-1]
             for k in ('nod_pntr', 'type', 'se_key', 'reluc', 'reluc_2',
                       'mag_1', 'mag_2', 'loss_dens')]
        logger.debug('Elements: %d', len(self.ELEM_ISA_ELEM_REC_EL_TYP))
             
        grp = ds.groups['element_nodes']
        (self.ELE_NOD_ISA_ND_KEY,
         self.ELE_NOD_ISA_NXT_ND_PNTR) = [
             grp.variables[k][:-1]
             for k in ('nd_key', 'nxt_nd_pntr')]
                        
        grp = ds.groups['superelements']
        (self.SUPEL_ISA_SE_NDCHN_PNTR,
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
         self.SUPEL_ISA_SUPEL_REC_SE_CURD_IM) = [
            grp.variables[k][:]
            for k in ('ndch_pntr', 'el_pntr', 'color', 'mcv_type',
                      'cond_type', 'vel_sys', 'sr_key', 'velo_1',
                      'velo_2', 'conduc', 'length', 'curd_re', 'curd_im')]
 
        grp = ds.groups['superelement_nodechains']
        (self.SE_NDCHN_ISA_NC_KEY,
         self.SE_NDCHN_ISA_NXT_NC_PNTR) = [
            grp.variables[k][:]
            for k in ('nc_key', 'nxt_nc_pntr')]
        
        grp = ds.groups['superelement_elements']
        (self.SE_EL_ISA_EL_KEY,
         self.SE_EL_ISA_NXT_EL_PNTR) = [
            grp.variables[k][:]
            for k in ('el_key', 'nxt_el_pntr')]

        grp = ds.groups['subregions']
        (self.SR_ISA_SR_SE_PNTR,
         self.SR_ISA_SR_REC_SR_TYP,
         self.SR_ISA_SR_REC_SR_COL,
         names,
         self.SR_ISA_SR_REC_SR_CUR_DIR,
         self.SR_ISA_SR_REC_SR_WB_KEY,
         self.SR_ISA_SR_REC_SR_NTURNS,
         self.SR_ISA_SR_REC_SR_SV_PNTR,
         self.SR_ISA_SR_REC_SR_ARRAY,
         self.SR_ISA_SR_REC_SR_GCUR_RE,
         self.SR_ISA_SR_REC_SR_GCUR_IM,
         self.SR_ISA_SR_REC_SR_VOLT_RE,
         self.SR_ISA_SR_REC_SR_VOLT_IM) = [
            grp.variables[k][:]
            for k in ('se_pntr', 'type', 'color', 'name', 'cur_dir',
                      'wb_key', 'nturns', 'sv_pntr', 'array',
                      'gcur_re', 'gcur_im', 'volt_re', 'volt_im')]
        self.SR_ISA_SR_REC_SR_NAME = [''.join(str(n, encoding='utf-8'))
                                      for n in names]
        grp = ds.groups['subregion_superelements']
        (self.SR_SE_ISA_SE_KEY,
         self.SR_SE_ISA_NXT_SE_PNTR) = [
            grp.variables[k][:]
            for k in ('se_key', 'nxt_se_pntr')]
        try:
            grp = ds.groups['windings']
            (self.WB_ISA_WB_SR_PNTR,
             self.WB_ISA_WB_REC_WB_COL,
             names,
             self.WB_ISA_WB_REC_WB_TURN,
             self.WB_ISA_WB_REC_WB_SR_NUM,
             self.WB_ISA_WB_REC_WB_WND_KEY,
             self.WB_ISA_WB_REC_WB_GCUR_RE,
             self.WB_ISA_WB_REC_WB_GCUR_IM,
             self.WB_ISA_WB_REC_WB_VOLT_RE,
             self.WB_ISA_WB_REC_WB_VOLT_IM,
             self.WB_ISA_WB_REC_WB_IMPDZ_RE,
             self.WB_ISA_WB_REC_WB_IMPDZ_IM) = [
                 grp.variables[k][:]
                 for k in ('sr_pntr', 'color', 'name', 'turn', 'sr_num',
                           'wnd_key', 'gcur_re', 'gcur_im',
                           'volt_re', 'volt_im', 'impdz_re', 'impdz_im')]
            self.WB_ISA_WB_REC_WB_UNIT_RES = self.WB_ISA_WB_REC_WB_TURN
            self.WB_ISA_WB_REC_WB_NAME = [''.join(str(n, encoding='utf-8'))
                                          for n in names]
        except:
            pass
        try:
            grp = ds.groups['winding_subregions']
            (self.WB_SR_ISA_SR_KEY,
             self.WB_SR_ISA_NXT_SR_PNTR) = [
                 grp.variables[k][:]
                 for k in ('sr_key', 'nxt_sr_pntr')]
        except:
            pass
        try:
            grp = ds.groups['machine']
            self.FC_RADIUS = float(grp.variables['fc_radius'].getValue().data)
            self.pole_pairs = int(grp.variables['pole_pairs'].getValue().data)
            self.poles_sim = int(grp.variables['poles_sim'].getValue().data)
            self.slots = int(grp.variables['num_slots'].getValue().data)
            self.arm_length = int(grp.variables['arm_length'].getValue().data)
        except:
            pass
        self.MAGN_TEMPERATURE = 20
        self.BR_TEMP_COEF = 0
        try:
            grp = ds.groups['magnet']
            self.MAGN_TEMPERATURE = float(grp.variables['temperature'].getValue().data)
            self.BR_TEMP_COEF = float(grp.variables['br_temp_coef'].getValue().data)
        except:
            pass
        try:
            grp = ds.groups['el_induction']
            (self.curr_loss,
             self.beta_loss,
             self.pos_el_fe_induction,
             self.el_fe_induction_1,
             self.el_fe_induction_2,
             self.eddy_cu_vpot) = [grp.variables[k][:]
                                   for k in ('cur', 'beta', 'position',
                                             'fe_induction_1',
                                             'fe_induction_2',
                                             'eddy_cu_vpot')]
            logger.debug('el_fe_induction %d', len(self.pos_el_fe_induction))
        except KeyError:
            self.pos_el_fe_induction = []
            self.el_fe_induction_1 = []
            self.el_fe_induction_2 = []
            self.eddy_cu_vpot = []

def read(filename):
    """
    Read nc file and return NcModel object.

    Arguments:
        filename: name of nc file to be read
    """
    import os
    ext = os.path.splitext(filename)[-1]
    if not ext:
        filename += '.nc'
 
    return isa7.Isa7(Reader(filename))


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s')
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    model = read(filename)
