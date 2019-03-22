# -*- coding: utf-8 -*-
"""
    femagtools.dxfsl.fslreader
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    renders geometry to fsl


"""
import numpy as np
from .shape import Shape
import logging
from .. import __version__
logger = logging.getLogger(__name__)


def gcd(a, b):
    """calc greatest common divisor"""
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """least common multiple"""
    return a*b//gcd(a, b)

    
def agndst2(da1, da2, Q, p):
    """return agndst from lcm"""
    r = (da1 + da2)/4
    ag = abs(da1 - da2)/6
    n = 12
    a = [r*4*np.pi/p/(i*lcm(Q, p))
         for i in range(1, n)]
    i = np.argmin(np.abs(np.array([ag]*(n-1)) - a))
    if i > 0:
        return a[i-1]
    return a[0]


def agndst(da1, da2, Q, p):
    """ build agndst from set of useful node angles:  
        4°, 2°, 1.5°, 1°, 0.75°, 0.5°, 0.25, 0.1, 0.05° """
    dagset = [np.pi/45, np.pi/90, np.pi/120, np.pi/180,
              np.pi/240, np.pi/360, np.pi/720, np.pi/1800, np.pi/3600]
    r = (da1 + da2)/4
    ag = abs(da1 - da2)/6
    i = max(np.argmin(np.abs(np.array(dagset) - np.arctan2(ag, r))), 1)
    return dagset[i-1]*r


class FslRenderer(object):
    """a model that can created by FSL"""

    def __init__(self, name):
        self.model = name
        self.mirror_axis = None
        self.fm_nlin = None
        self.shaft = None

    def mirror_nodechains(self, p0, p1):
        self.mirror_axis = np.array((p0, p1)).ravel().tolist()

    def material(self, p0):
        self.fm_nlin = p0

    def circle(self, center, radius, color='blue'):
        num = int(2*np.pi*radius)
        if num < 8:
            num = 8
        circle = [u'cx, cy = {}, {}'.format(center[0],
                                            center[1]),
                  u'nc_circle_m(cx,cy,{}, {})'.format(radius, num),
                  u'create_mesh_se(cx, cy)\n']
        self.content += circle

    def arc(self, startangle, endangle, center, radius, color='blue'):
        num = 0
#        if self.nodedist > 0:
#            s = startangle
#            d = endangle - s
#            n = int(d/(2*np.pi))
#            if n < 0:
#                n -= 1
#            d -= n*2*np.pi
#            num = int(radius*d/self.nodedist + 1)
#            if num < 3 and radius*d > self.nodedist:
#                num = 3

        p1 = (center[0] + radius*np.cos(startangle),
              center[1] + radius*np.sin(startangle))
        p2 = (center[0] + radius*np.cos(endangle),
              center[1] + radius*np.sin(endangle))
        self.content.append(
            u"nc_circle_m({}, {}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1], num))

    def line(self, p1, p2, color='blue', e=None):
        num = 0
        # if self.nodedist > 0:
        #     l = la.norm(np.asarray(p1)-p2)
        #     num = int(l/self.nodedist + 1)
        if e is not None and e.has_attribute('auxline'):
            self.content.append(
                u"nc_line({}, {}, {}, {}, {}) -- auxiliary".format(
                    p1[0], p1[1], p2[0], p2[1], num))
        else:
            self.content.append(
                u"nc_line({}, {}, {}, {}, {})".format(
                    p1[0], p1[1], p2[0], p2[1], num))

    def sorted_elements(self, geom, inner=False):
        if inner:
            # Abstand von airgap Richtung Nullpunkt
            el_sorted = [(geom.max_radius -
                          e.minmax_from_center((0.0, 0.0))[1], e)
                         for e in geom.elements(Shape)]
        else:
            # Abstand von airgap Richtung Aussen
            el_sorted = [(e.minmax_from_center((0.0, 0.0))[0] -
                          geom.min_radius, e)
                         for e in geom.elements(Shape)]

        el_sorted.sort()
        return el_sorted

    def render(self, machine, inner=False, outer=False):
        '''create fsl statements with nodechains'''
        machine.set_alfa_and_corners()
        geom = machine.geom
        self.content = []

        ndt_list = [(0.25, 1.5), (0.5, 2), (0.75, 3.0), (1.1, 3.0)]
        dist = geom.max_radius - geom.min_radius
        el_sorted = self.sorted_elements(geom, inner)

        x = 0
        for d, e in el_sorted:
            d_percent = d / dist
            if ndt_list[x][0] < d_percent:
                self.content.append(u'\nndt({}*agndst)\n'.
                                    format(ndt_list[x][1]))
                while len(ndt_list) and ndt_list[x][0] < d_percent:
                    x += 1
#            self.content.append(u'-- d={} / dist={} == {}'.
#                                format(d, dist, d_percent))
            e.render(self)

        self.content.append(u'\n')

        parts = int(machine.get_symmetry_part())
        self.content += [u'-- parts      = {}'.format(parts),
                         u'-- min_radius = {}'.format(geom.min_radius),
                         u'-- max_radius = {}'.format(geom.max_radius),
                         u'-- min_corner = {}, {}'.format(
                             geom.start_min_corner(0),
                             geom.start_min_corner(1)),
                         u'-- max_corner = {}, {}'.format(
                             geom.start_max_corner(0),
                             geom.start_max_corner(1)),
                         u'\n']
        if inner:
            self.content += [u'inner_da_start = {}'
                             .format(geom.dist_start_max_corner()),
                             u'inner_da_end = {}'
                             .format(geom.dist_end_max_corner())]
        if outer:
            self.content += [u'-- create air layer outside',
                             u'x0, y0 = {}, {}'.format(
                                 geom.start_max_corner(0),
                                 geom.start_max_corner(1)),
                             u'hair = 1.0',
                             u'r1 = {} + hair'.format(geom.max_radius),
                             u'r, phi = c2pr(x0, y0)',
                             u'x1, y1 = pr2c(r1, phi)']
            if geom.is_mirrored():
                self.content += [
                    u'x2, y2 = pr2c(r1, math.pi/m.tot_num_slot+phi/2)',
                    u'x3, y3 = pr2c(r, math.pi/m.tot_num_slot+phi/2)']
            else:
                self.content += [
                    u'x2, y2 = pr2c(r1, 2*math.pi/m.tot_num_slot+phi)',
                    u'x3, y3 = pr2c(r, 2*math.pi/m.tot_num_slot+phi)']
                
            self.content += [u'nc_line(x0, y0, x1, y1, 0)',
                             u'nc_circle_m(x1, y1, x2, y2, 0.0, 0.0, 0)',
                             u'nc_line(x2, y2, x3, y3, 0)']
            if geom.is_mirrored():
                self.content.append(
                    u'x0, y0 = pr2c(r1 - hair/2, math.pi/m.tot_num_slot/2+phi/4)')
            else:    
                self.content.append(
                    u'x0, y0 = pr2c(r1 - hair/2, math.pi/m.tot_num_slot+phi/2)')
                
            self.content += [
                u'create_mesh_se(x0, y0)',
                u'\n',
                u'outer_da_start = {}'.format(
                    geom.dist_start_min_corner()),
                u'outer_da_end = {}'.format(
                    geom.dist_end_min_corner())]

        self.content.append(u'\n')
        self.content.append(u'x0_iron_shaft, y0_iron_shaft = 0.0, 0.0')
        self.content.append(u'x0_iron_yoke, y0_iron_yoke = 0.0, 0.0\n')

        subregions = {}
        num_windings = 0
        num_magnets = 0
        for area in geom.list_of_areas():
            if area.number_of_elements() > 1:
                p = area.get_point_inside(geom)
                if p:
                    self.content.append(u"x0, y0 = {}, {}".format(p[0], p[1]))
                    self.content.append(u"point(x0, y0, red, 4)")
                    self.content.append(u"create_mesh_se(x0, y0)")

                if area.is_winding():
                    if area.type not in subregions:
                        subregions[area.type] = 1
                    num_windings += 1
                    self.content.append(u'm.xcoil_{}, m.ycoil_{} = x0, y0'.
                                        format(num_windings, num_windings))

                elif area.is_magnet():
                    if area.type not in subregions:
                        subregions[area.type] = 1
                    num_magnets += 1
                    self.content.append(u'xmag[{}], ymag[{}] = x0, y0'.
                                        format(num_magnets, num_magnets))
                    self.content.append(u'mag_orient[{}] = {}'.
                                        format(num_magnets, area.phi))
                elif area.type > 0:
                    if area.type in subregions:
                        self.content.append(
                            u'add_to_subreg(x0, y0, "{}")'.
                            format(area.name()))
                    else:
                        subregions[area.type] = 1
                        self.content.append(
                            u'def_new_subreg(x0, y0, "{}", {})'.
                            format(area.name(), area.color()))
                    if area.is_stator_iron_yoke():
                        self.content.append(
                            u'x0_iron_yoke, y0_iron_yoke = x0, y0')
                    if area.is_stator_iron_shaft():
                        self.content.append(
                            u'x0_iron_shaft, y0_iron_shaft = x0, y0')
                    if area.is_rotor_iron():
                        self.content.append(
                            u'x0_iron_yoke, y0_iron_yoke = x0, y0')

                self.content.append(u"\n")

        txt = [u"if x0_iron_yoke > 0.0 then",
               u"  if mcvkey_yoke ~= 'dummy' then",
               u'    def_mat_fm_nlin(x0_iron_yoke, y0_iron_yoke, blue, mcvkey_yoke, 100)',
               u'  else',
               u'    def_mat_fm(x0_iron_yoke, y0_iron_yoke, ur, 100)',
               u'  end',
               u'end\n']
        self.content.append(u'\n'.join(txt))

        txt = [u"if x0_iron_shaft > 0.0 then",
               u"  if mcvkey_shaft ~= 'dummy' then",
               u'    def_mat_fm_nlin(x0_iron_shaft, y0_iron_shaft, blue, mcvkey_shaft, 100)',
               u'  else',
               u'    def_mat_fm(x0_iron_shaft, y0_iron_shaft, ur, 100)',
               u'  end',
               u'end\n']
        self.content.append(u'\n'.join(txt))

        if num_windings > 0:
            if geom.is_mirrored():
                self.content += [
                    u'    r, phi = c2pr(m.xcoil_1, m.ycoil_1)',
                    u'    m.xcoil_2, m.ycoil_2 = pr2c(r, {}*2.0 - phi)'.format(geom.alfa)]
            self.content.append(u'm.wdg_location  = 1.0 -- stator\n')

        if num_magnets > 0:
            self.content.append(u'mag_exists   = {}'.
                                format(num_magnets))
            if geom.is_mirrored():
                self.content.append(u'mag_mirrored = true')
            else:
                self.content.append(u'mag_mirrored = false')
            self.content.append(u'mag_alpha    = {}\n'
                                .format(geom.get_alfa()))

        if geom.is_mirrored():
            self.content.append(u'-- mirror')
            self.content.append(u'mirror_nodechains({}, {}, {}, {})\n'.format(
                geom.mirror_corners[1][0],   # max x1
                geom.mirror_corners[1][1],   # max y1
                geom.mirror_corners[0][0],   # min x2
                geom.mirror_corners[0][1]))  # min y2

        # angle after mirroring
        self.content.append(u'alfa = {}\n'.format(geom.get_alfa()))

        self.content += [u'-- rotate',
                         u'x1, y1 = {}, {}'.format(
                             geom.start_corners[0][0],
                             geom.start_corners[0][1])] # min xy1
        if outer:
            self.content.append(u'x2, y2 = pr2c(r1, phi)')
        else:
            self.content.append(u'x2, y2 = {}, {}'.format(
                             geom.start_corners[1][0],
                             geom.start_corners[1][1])) # max xy1
            
        if geom.is_mirrored():
            self.content.append(u'x3, y3 = pr2c(x2, alfa)')
            self.content.append(u'x4, y4 = pr2c(x1, alfa)')
        else:
            self.content += [u'x3, y3 = pr2c(r1, 2*math.pi/m.tot_num_slot+phi)',
                             u'x4, y4 = {}, {}'.format(
                                 geom.end_corners[0][0],
                                 geom.end_corners[0][1])]  # min xy4

        self.content.append(u'if {} > 1 then'.format(geom.num_variable()))
        if geom.corners_dont_match():
            txt = [u'  -- Warning: corners dont match',
                   u'  mirror_nodechains(x3, y3, x4, y4)',
                   u'  if {} > 2 then'.format(geom.num_variable()),
                   u'    my_alfa = alfa * 2',
                   u'    x3, y3 = pr2c(x2, my_alfa)',
                   u'    x4, y4 = pr2c(x1, my_alfa)',
                   u'    n = {} / 2'.format(geom.num_variable()),
                   u'    rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,n-1)',
                   u'  end']
            self.content.append(u'\n'.join(txt))
        else:
            self.content.append(
                u'  rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,{}-1)'
                .format(geom.num_variable()))
        self.content.append(u'end')

        self.content.append(u'alfa = {} * alfa'.format(geom.num_variable()))

        if self.fm_nlin:
            self.content.append(u'\nx0, y0 = {}, {}'. format(
                self.fm_nlin[0], self.fm_nlin[1]))

            mat = [u"if fm_nlin_mcvfile ~= 'dummy' then",
                   u"  if fm_nlin_mcvfile == 'air' then",
                   u"    def_mat_fm(x0,y0, 1.0, fm_nlin_rlen)",
                   u"  else",
                   u"    def_mat_fm_nlin(x0,y0, fm_nlin_colour, fm_nlin_mcvfile, fm_nlin_rlen)",
                   u"  end",
                   u"else",
                   u"  def_mat_fm(x0,y0, 1000.0, fm_nlin_rlen)",
                   u"end"]
            self.content.append(u'\n'.join(mat))

        if self.shaft:
            mat = [u'\nif shaft_mat==1 then',
                   u'  def_mat_fm_nlin(0.1,0.1,lightgrey,fm_nlin_mcvfile_shft,fm_nlin_rlen)',
                   u'end']
            self.content.append(u'\n'.join(mat))

        if num_magnets:
            self.content += [
                u'-- pm magnets',
                u'if mag_exists > 0 then',
                u'  m.remanenc = 1.15',
                u'  m.relperm  = 1.05',
                u'  for i=0, m.npols_gen-1 do',
                u'    for n=1, mag_exists do',
                u'      r, p = c2pr(xmag[n], ymag[n])',
                u'      phi = i*mag_alpha+p',
                u'      x0, y0 = pr2c(r, phi)',
                u'      phi_orient = i*mag_alpha+mag_orient[n]',
                u'      orient = phi_orient*180/math.pi',
                u'      if ( i % 2 == 0 ) then',
                u'        def_mat_pm(x0, y0, red, m.remanenc, m.relperm,',
                u'                   orient-180, m.parallel, 100)',
                u'      else',
                u'        def_mat_pm(x0, y0, green, m.remanenc, m.relperm,',
                u'                   orient, m.parallel, 100)',
                u'      end',
                u'      if mag_mirrored then',
                u'        phi = (i+1)*mag_alpha-p',
                u'        x0, y0 = pr2c(r, phi)',
                u'        phi_orient = (i+1)*mag_alpha-mag_orient[n]',
                u'        orient = phi_orient*180/math.pi',
                u'        if ( i % 2 == 0 ) then',
                u'          def_mat_pm(x0, y0, red, m.remanenc, m.relperm,',
                u'                     orient-180, m.parallel, 100)',
                u'        else',
                u'          def_mat_pm(x0, y0, green, m.remanenc, m.relperm,',
                u'                     orient, m.parallel, 100)',
                u'        end',
                u'      end',
                u'    end',
                u'  end',
                u'end']

        return self.content

    def render_main(self,
                    motor,
                    m_inner, m_outer,
                    inner, outer,
                    params):
        '''create fsl statements for main file'''
        if not (m_inner and m_outer):
            logger.warning("ERROR: Rotor or Stator missing")
            return []

        num_layers = min(m_inner.num_of_layers() +
                         m_outer.num_of_layers(),
                         2)
        return [
            u'-- generated from DXF by femagtools {}'.format(__version__),
            u'exit_on_error = false',
            u'exit_on_end = false',
            u'verbosity = 2\n',
            u'new_model_force("{}","Test")'.format(self.model),
            u'global_unit(mm)',
            u'pickdist(0.001)',
            u'cosys(polar)\n',
            u'xmag = {}',
            u'ymag = {}',
            u'mag_orient = {}',
            u'mag_exists = 0',
            u'dy1 = {}'.format(params.get('dy1', 0.0)),
            u'da1 = {}'.format(params.get('da1', 0.0)),
            u'dy2 = {}'.format(params.get('dy2', 0.0)),
            u'da2 = {}'.format(params.get('da2', 0.0)),
            u'ag  = (da1 - da2)/2\n',
            u'm.tot_num_slot   = {}'.format(params.get('tot_num_slot', 0)),
            u'm.num_sl_gen     = {}'.format(params.get('num_sl_gen', 0)),
            u'm.num_poles      = {}'.format(params.get('num_poles', 0)),
            u'm.num_pol_pair   = m.num_poles/2',
            u'm.num_slots      = m.num_sl_gen',
            u'm.npols_gen      = m.num_poles * m.num_sl_gen / m.tot_num_slot',
            u'm.tot_num_sl     = m.tot_num_slot',
            u'm.fc_radius      = (da1+da2)/4',
            u'm.fc_radius1     = m.fc_radius',
            u'm.arm_length     = 1.0',
            u'pre_models("basic_modpar")\n',
            u'm.airgap         = 2*ag/3',
            u'm.nodedist       = 1.0',
            u'agndst           = {}'.format(params.get('agndst', 0.1)),
            u"mcvkey_yoke = 'dummy'",
            u"mcvkey_shaft = 'dummy'",
            u"ur = 1000.0",
            u"ndt(agndst)"] + outer + [
                "mcvkey_yoke = 'dummy'",
                u"mcvkey_shaft = 'dummy'",
                u"ndt(agndst)"] + inner + [
                    u'-- airgap',
                    u'ndt(agndst)',
                    u'r1 = da2/2 + ag/3',
                    u'x1, y1 = pr2c(r1, alfa)',
                    u'n = r1*alfa/agndst + 1',
                    u'nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)\n',
                    u'r2 = da2/2 + 2*ag/3',
                    u'x2, y2 = pr2c(r2, alfa)',
                    u'nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)\n',
                    u'if inner_da_start == nil then',
                    u'  inner_da_start = da2/2',
                    u'end',
                    u'x1, y1 = pr2c(inner_da_start, 0.0)',
                    u'nc_line(x1, y1, r1, 0.0, 0.0)\n',
                    u'if outer_da_start == nil then',
                    u'  outer_da_start = da1/2',
                    u'end',
                    u'x2, y2 = pr2c(outer_da_start, 0.0)',
                    u'nc_line(r2, 0.0, x2, y2, 0.0)\n',
                    u'if m.tot_num_slot > m.num_sl_gen then',
                    u'  x3, y3 = pr2c(inner_da_end, alfa)',
                    u'  x4, y4 = pr2c(r1, alfa)',
                    u'  nc_line(x3, y3, x4, y4, 0, 0)\n',
                    u'  x3, y3 = pr2c(outer_da_end, alfa)',
                    u'  x4, y4 = pr2c(r2, alfa)',
                    u'  nc_line(x3, y3, x4, y4, 0, 0)',
                    u'end\n',
                    u'x0, y0 = pr2c(r1-ag/6, alfa/2)',
                    u'create_mesh_se(x0, y0)',
                    u'x0, y0 = pr2c(r2+ag/6, alfa/2)',
                    u'create_mesh_se(x0, y0)\n',
                    u'connect_models()\n',
                    u'-- Gen_winding',
                    u'if m.xcoil_1 ~= nil then',
                    u'  m.num_phases      = 3',
                    u'  m.num_layers      = {}'.format(num_layers),
                    u'  m.num_wires       = 1',
                    u'  m.coil_span       = 1',
                    u'  m.current         = 0.0',
                    u'  m.mat_type        = 1.0 -- rotating',
                    u'  m.wind_type       = 1.0 -- winding & current',
                    u'  m.win_asym        = 1.0 -- sym',
                    u'  m.wdg_location    = 1.0 -- stator',
                    u'  m.curr_inp        = 0.0 -- const',
                    u'  m.dq_offset       = 0',
                    u'  pre_models("Gen_winding")',
                    u'  pre_models("gen_pocfile")',
                    u'end\n',
                    u'save_model(cont)\n']
