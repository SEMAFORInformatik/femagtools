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


def agndst(da1, da2, Q, p, nodedist=1):
    """ build agndst from set of useful number of nodes"""
    num_nodes = [30, 48, 60, 96, 120, 144, 180, 240, 288, 336, 360,
                 432, 480]
    r = (da1 + da2)/4
    dagset = [2*np.pi/p/i for i in num_nodes]
    ag = abs(da1 - da2)/6
    i = max(np.argmin(np.abs(np.array(dagset) - np.arctan2(ag, r))), 1)
    nd = min(round(nodedist), i)
    try:
        logger.info("Num nodes/p %d Num nodes/slot %g nodedist %g",
                    num_nodes[i-1], p*num_nodes[i-1]/Q, nodedist)
        if nodedist > 1:
            return dagset[i-nd]*r
        elif nodedist < 1 or i == 0:
            return dagset[i]*r
    except IndexError:
        pass
    return dagset[i-1]*r


class FslRenderer(object):
    """a model that can created by FSL"""

    def __init__(self, name):
        self.model = name
        self.mirror_axis = None
        self.fm_nlin = None
        self.shaft = None
        self.agndst = 1

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
        #d = (endangle - startangle) % 2*np.pi
        #num = int(radius*d/self.agndst + 1)
        # if num < 3 and d > np.pi/6:
        #    num = 3
        # else:
        #    num = 0

        p1 = (center[0] + radius*np.cos(startangle),
              center[1] + radius*np.sin(startangle))
        p2 = (center[0] + radius*np.cos(endangle),
              center[1] + radius*np.sin(endangle))
        self.content.append(
            u"nc_circle_m({}, {}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1], num))

    def ellipse(self, center, width, height,
                rtheta, start_param, end_param, color='blue'):
        theta = np.arange(start_param, end_param, 0.2)
        xbegin = 0.5 * width * np.cos(theta)
        ybegin = 0.5 * height * np.sin(theta)
        theta = np.arange(end_param, end_param+0.1, 0.3)
        xend = 0.5 * width * np.cos(theta)
        yend = 0.5 * height * np.sin(theta)

        x = np.concatenate((xbegin, xend), axis=0)
        y = np.concatenate((ybegin, yend), axis=0)

        R = np.array([
            [np.cos(rtheta), -np.sin(rtheta)],
            [np.sin(rtheta),  np.cos(rtheta)]])
        x, y = np.dot(R, np.array([x, y]))
        x += center[0]
        y += center[1]
        nodes = list(zip(x, y))
        n0 = nodes[0]
        for n1 in nodes[1:]:
            self.content.append(
                u"nc_line({}, {}, {}, {}, {}) -- ellipse".format(
                    n0[0], n0[1], n1[0], n1[1], 0))
            n0 = n1

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
            r = geom.max_radius
        else:
            r = geom.min_radius

        # return sorted list of distances from airgap and elements
        return sorted([(abs(r - np.linalg.norm(e.center_of_connection())), e)
                       for e in geom.elements(Shape)])

    def render(self, machine, inner=False, outer=False):
        '''create fsl statements with nodechains'''
        machine.set_alfa_and_corners()
        geom = machine.geom
        geom.split_lines_longer_than(geom.max_radius/4)
        self.content = []

        ndt_list = [(0.2, 1.5), (0.45, 2), (0.7, 3.0), (1.1, 3.0)]
        dist = geom.max_radius - geom.min_radius
        el_sorted = self.sorted_elements(geom, inner)

        n = 0
        for d, e in el_sorted:
            d_percent = d / dist
            if ndt_list[n][0] < d_percent:
                self.agndst = ndt_list[n][1] * self.agndst
                self.content.append(u'\nndt({}*agndst)\n'.
                                    format(ndt_list[n][1]))
                while len(ndt_list) and ndt_list[n][0] < d_percent:
                    n += 1
            e.render(self)

        self.content.append(u'\n')

        parts = int(machine.get_symmetry_part())
        self.content += [u'-- parts      = {}'.format(parts)]
        if geom.is_stator():
            if machine.get_num_slots() > 0:
                self.content += [
                    u'-- num_slots  = {}'.format(
                        machine.get_num_slots())
                ]
        self.content += [u'-- min_radius = {}'.format(geom.min_radius),
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
            self.content += [
                u'-- create air layer outside',
                u'x0, y0 = {}, {}'.format(
                    geom.start_max_corner(0),
                    geom.start_max_corner(1)),
                u'hair = 1.0',
                u'parts = {}'.format(machine.get_num_parts()),
                u'r1 = {} + hair'.format(geom.max_radius),
                u'r, phi = c2pr(x0, y0)',
                u'x1, y1 = pr2c(r1, phi)',
                u'x2, y2 = pr2c(r1, 2*math.pi/parts)',
                u'x3, y3 = pr2c(r, 2*math.pi/parts)',
                u'nc_line(x0, y0, x1, y1, 0)',
                u'nc_circle_m(x1, y1, x2, y2, 0.0, 0.0, 0)',
                u'nc_line(x2, y2, x3, y3, 0)',
                u'x0, y0 = pr2c(r1 - hair/2, math.pi/parts)',
                u'create_mesh_se(x0, y0)',
                u'\n',
                u'outer_da_start = {}'.format(
                    geom.dist_start_min_corner()),
                u'outer_da_end = {}'.format(
                    geom.dist_end_min_corner())
            ]

        self.content += [u'\n',
                         u'xmag = {}',
                         u'ymag = {}',
                         u'mag_orient = {}',
                         u'mag_exists = 0',
                         u'if mcvkey_yoke == nil then',
                         u'  mcvkey_yoke = "dummy"',
                         u'end',
                         u'x0_iron_tooth, y0_iron_tooth = 0.0, 0.0',
                         u'x0_iron_yoke, y0_iron_yoke = 0.0, 0.0',
                         u'x0_shaft, y0_shaft = 0.0, 0.0',
                         u'\n']

        subregions = {}
        num_windings = 0
        num_magnets = 0
        for area in geom.list_of_areas():
            if area.number_of_elements() > 1:
                p = area.get_point_inside(geom)
                if p:
                    self.content.append(u"x0, y0 = {}, {}".format(p[0], p[1]))
                    # self.content.append(u"point(x0, y0, red, 4)")  # for debugging
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
                            u'def_new_subreg(x0, y0, "{}", "{}")'.
                            format(area.name(), area.color()))
                    if area.is_stator_iron_yoke():
                        self.content.append(
                            u'x0_iron_yoke, y0_iron_yoke = x0, y0')
                    elif area.is_stator_iron_tooth():
                        self.content.append(
                            u'x0_iron_tooth, y0_iron_tooth = x0, y0')
                    elif area.is_rotor_iron():
                        self.content.append(
                            u'x0_iron_yoke, y0_iron_yoke = x0, y0')
                    elif area.is_shaft():
                        self.content.append(
                            u'x0_shaft, y0_shaft = x0, y0')

                self.content.append(u"\n")

        txt = [u"if x0_iron_yoke > 0.0 then",
               u"  if mcvkey_yoke ~= 'dummy' then",
               u'    def_mat_fm_nlin(x0_iron_yoke, y0_iron_yoke, "blue", mcvkey_yoke, 100)',
               u'  else',
               u'    def_mat_fm(x0_iron_yoke, y0_iron_yoke, ur, 100)',
               u'  end',
               u'end\n']
        self.content.append(u'\n'.join(txt))

        txt = [u"if x0_iron_tooth > 0.0 then",
               u"  if(x0_iron_yoke == 0 and mcvkey_yoke ~= 'dummy') then",
               u"    def_mat_fm_nlin(x0_iron_tooth, y0_iron_tooth, 'blue', mcvkey_yoke, 100)",
               u"  else",
               u"    if (mcvkey_teeth ~= 'dummy' and mcvkey_teeth ~= nil) then",
               u"      def_mat_fm_nlin(x0_iron_tooth, y0_iron_tooth, 'blue', mcvkey_teeth, 100)",
               u"    else",
               u"      def_mat_fm(x0_iron_tooth, y0_iron_tooth, ur, 100)",
               u"    end",
               u"  end",
               u'end\n']
        self.content.append(u'\n'.join(txt))

        txt = [u"if x0_shaft > 0.0 then",
               u"  if mcvkey_shaft ~= 'dummy' then",
               u'    def_mat_fm_nlin(x0_shaft, y0_shaft, "lightgrey", mcvkey_shaft, 100)',
               u'  else',
               u'    def_mat_fm(x0_shaft, y0_shaft, ur, 100)',
               u'  end',
               u'end\n']
        self.content.append(u'\n'.join(txt))

        if num_windings > 0:
            if geom.is_mirrored():
                self.content += [
                    u'r, phi = c2pr(m.xcoil_1, m.ycoil_1)',
                    u'm.xcoil_2, m.ycoil_2 = pr2c(r, {}*2.0 - phi)'.format(geom.alfa)]
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

        num_parts = machine.get_num_parts()
        if num_parts == parts:
            self.content.append(u'parts_gen = {}'.format(geom.num_variable()))
        else:
            self.content.append(u'parts_gen = {}/2'.format(geom.num_variable()))

        # angle after mirroring
        self.content.append(u'alfa = {}\n'.format(geom.get_alfa()))

        self.content += [u'-- rotate',
                         u'x1, y1 = {}, {}'.format(
                             geom.start_corners[0][0],
                             geom.start_corners[0][1])]  # min xy1
        if outer:
            self.content.append(u'x2, y2 = pr2c(r1, 0.0)')
        else:
            self.content.append(u'x2, y2 = {}, {}'.format(
                geom.start_corners[1][0],
                geom.start_corners[1][1]))  # max xy1

        if geom.is_mirrored():
            self.content.append(u'x3, y3 = pr2c(x2, alfa)')
            self.content.append(u'x4, y4 = pr2c(x1, alfa)')
        else:
            self.content += [u'x3, y3 = pr2c(r1, 2*math.pi/m.tot_num_slot+phi)',
                             u'x4, y4 = {}, {}'.format(
                                 geom.end_corners[0][0],
                                 geom.end_corners[0][1])]  # min xy4


        self.content.append(u'if parts_gen > 1 then')
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
                u'  rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,parts_gen-1)')
        self.content.append(u'end')

        self.content.append(u'alfa = parts_gen * alfa\n')

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
                   u'  def_mat_fm_nlin(0.1,0.1,"lightgrey",fm_nlin_mcvfile_shft,fm_nlin_rlen)',
                   u'end']
            self.content.append(u'\n'.join(mat))

        if num_magnets:
            self.content += [
                u'-- pm magnets',
                u'if mag_exists > 0 then',
                u'  if( m.remanenc == nil ) then',
                u'     m.remanenc = 1.15',
                u'     m.relperm  = 1.05',
                u'  end',
                u'  if( m.orient == nil ) then',
                u'    m.orient = m.parallel',
                u'  end',
                u'  if( m.magncond == nil ) then',
                u'    m.magncond = 625000.0',
                u'  end',
                u'  if( m.rlen == nil ) then',
                u'    m.rlen = 100',
                u'  end',
                u'  for i=0, m.npols_gen-1 do',
                u'    for n=1, mag_exists do',
                u'      r, p = c2pr(xmag[n], ymag[n])',
                u'      x0, y0 = pr2c(r, i*mag_alpha+p)',
                u'      phi = (i*mag_alpha+mag_orient[n])*180/math.pi',
                u'      if ( i % 2 == 0 ) then',
                u'        phi = phi - 180',
                u'        color = "red"',
                u'      else',
                u'        color = "green"',
                u'      end',
                u'      if(m.mcvkey_magnet == nil) then',
                u'        def_mat_pm(x0, y0, color, m.remanenc, m.relperm,',
                u'                   phi, m.orient, m.magncond, m.rlen)',
                u'      else',
                u'        def_mat_pm_nlin(x0, y0, color, m.mcvkey_magnet,',
                u'                 phi, m.orient, m.magncond, m.rlen)',
                u'      end',
                u'      if mag_mirrored then',
                u'        x0, y0 = pr2c(r, (i+1)*mag_alpha-p)',
                u'        phi = ((i+1)*mag_alpha-mag_orient[n])*180/math.pi',
                u'        if ( i % 2 == 0 ) then',
                u'          phi = phi - 180',
                u'        end',
                u'        if(m.mcvkey_magnet == nil) then',
                u'          def_mat_pm(x0, y0, color, m.remanenc, m.relperm,',
                u'                       phi, m.orient, m.magncond, m.rlen)',
                u'        else',
                u'          def_mat_pm_nlin(x0, y0, color, m.mcvkey_magnet,',
                u'                   phi, m.orient, m.magncond, m.rlen)',
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
        self.agndst = params.get('agndst', 0.1)
        return [
            u'-- generated from DXF by femagtools {}'.format(__version__),
            u'exit_on_error = false',
            u'exit_on_end = false',
            u'verbosity = 2\n',
            u'new_model_force("{}","Test")'.format(self.model),
            u'global_unit(mm)',
            u'pickdist(0.001)',
            u'cosys(polar)\n',
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
            u'agndst           = {}'.format(self.agndst),
            u"mcvkey_teeth = 'dummy'",
            u"mcvkey_yoke = 'dummy'",
            u"mcvkey_shaft = 'dummy'",
            u"ur = 1000.0",
            u"ndt(agndst)"] + outer + [
                u"mcvkey_teeth = 'dummy'",
                u"mcvkey_yoke = 'dummy'",
                u"mcvkey_shaft = 'dummy'",
                u"ndt(agndst)"] + inner + [
                    u'-- airgap',
                    u'ndt(agndst)',
                    u'r1 = m.fc_radius - ag/6',
                    u'x1, y1 = pr2c(r1, alfa)',
                    u'n = math.floor(r1*alfa/agndst + 1.5)',
                    u'nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)\n',
                    u'r2 = m.fc_radius + ag/6',
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
                    u'  m.coil_span       = {}'.format(
                        2*num_layers*params.get(
                            'tot_num_slot', 1)//params.get('num_poles', 1)//3),
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
