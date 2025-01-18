# -*- coding: utf-8 -*-
"""
    femagtools.dxfsl.fslreader
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    renders geometry to fsl


"""
import numpy as np
from .shape import Shape
from .area import TYPE_TOOTH, TYPE_YOKE
import logging
from .. import __version__
logger = logging.getLogger(__name__)


def agndst(da1, da2, Q, p, nodedist=1):
    """ build agndst from set of useful number of nodes
    args:
      da1: stator inner diameter
      da2: rotor outer diameter
      Q: number of stator slots
      p: number of poles
      nodedist: node distance factor
    """
    r = (da1 + da2)/4
    ag = abs(da1 - da2)/6
    lcm = np.lcm(Q, p)//p
    taup = 2*np.pi*r/p
    nmin, nmax = max(1, 18//lcm), 480//lcm
    num_nodes = [i*lcm for i in range(nmin, nmax)
                 if i*lcm % 6 == 0 and i*lcm*p//Q % 2 == 0]
    dagset = np.array([taup/i for i in num_nodes])
    i = np.argmin(np.abs(dagset - ag*nodedist))
    if i > 0 and dagset[i]*nodedist < ag:
        i -= 1
    logger.info("%d Num nodes/p %d Num nodes/slot %g nodedist %g",
                i, num_nodes[i], p*num_nodes[i]/Q, nodedist)
    return dagset[i]


class FslRenderer(object):
    """a model that can created by FSL"""

    def __init__(self, name, mtype):
        self.model = name
        self.mirror_axis = None
        self.fm_nlin = None
        self.shaft = None
        self.agndst = 1
        self.mtype = mtype

    def mirror_nodechains(self, p0, p1):
        self.mirror_axis = np.array((p0, p1)).ravel().tolist()

    def material(self, p0):
        self.fm_nlin = p0

    def circle(self, center, radius,
               color='blue', linestyle=None):
        num = int(2*np.pi*radius)
        if num < 8:
            num = 8
        circle = ['cx, cy = {}, {}'.format(center[0],
                                            center[1]),
                  'nc_circle_m(cx,cy,{}, {})'.format(radius, num),
                  'create_mesh_se(cx, cy)\n']
        self.content += circle

    def arc(self, startangle, endangle, center, radius,
            color='blue', linestyle=None):
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
            "nc_circle_m({}, {}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1], num))

    def ellipse(self, center, width, height,
                rtheta, start_param, end_param,
                color='blue', linestyle=None):
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
                "nc_line({}, {}, {}, {}, {}) -- ellipse".format(
                    n0[0], n0[1], n1[0], n1[1], 0))
            n0 = n1

    def line(self, p1, p2, e=None,
             color='blue', linestyle=None):
        num = 0
        # if self.nodedist > 0:
        #     l = la.norm(np.asarray(p1)-p2)
        #     num = int(l/self.nodedist + 1)
        if e is not None:
            if e.has_attribute('no_fsl'):
                logger.info("line with attr nofsl")
                return
            if e.has_attribute('auxline'):
                self.content.append(
                    "nc_line({}, {}, {}, {}, {}) -- auxiliary".format(
                        p1[0], p1[1], p2[0], p2[1], num))
                return

        self.content.append(
            "nc_line({}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1], num))

    def sorted_elements(self, geom, inner=False):
        if inner:
            r = geom.max_radius
        else:
            r = geom.min_radius

        # return sorted list of distances from airgap and elements
        return sorted([(abs(r - np.linalg.norm(e.center_of_connection())), e)
                       for e in geom.elements(Shape)])

    def render(self, machine, inner=False, outer=False, standalone=False):
        '''create fsl statements with nodechains'''
        machine.set_alfa_and_corners()
        geom = machine.geom
        split_len = (geom.max_radius - geom.min_radius) / 4
        geom.split_all_lines_longer_than(split_len)
        self.content = []

        if standalone:
            self.content += ['if (agndst == nil) then',
                             '  agndst = 0.5',
                             '  m.npols_gen = 2',
                             '  m.num_sl_gen = 2',
                             '  new_model_force("{}","Test")'.format(self.model),
                             'end']

        MAXDST=4.0
        NUMLEVELS=10
        NDT0=1.1
        # ndt list format [ (rdistx, ndtx) ...]
        # where
        #   - rdistx is rel dist from airgap (range 0 .. NUMLEVELS-1/NUMLEVELS)
        #   - ndtx nodedist (range NDT0 .. (NUMLEVELS-1/NUMLEVELS))*(MAXDST-1.1) + NDT0)
        ndt_list = [(1.1*nl/NUMLEVELS, nl/NUMLEVELS*(MAXDST-1.1)+NDT0)
                    for nl in range(NUMLEVELS+1)]

        dist = geom.max_radius - geom.min_radius
        el_sorted = self.sorted_elements(geom, inner)

        n = 0
        for d, e in el_sorted:
            d_percent = min(1.0, d / dist)
            if ndt_list[n][0] < d_percent:
                self.agndst = ndt_list[n][1] * self.agndst
                self.content.append('\nndt({}*agndst)\n'.
                                    format(round(ndt_list[n][1], 1)))
                while ndt_list[n][0] < d_percent:
                    n += 1
            e.render(self)

        self.content.append('\n')

        parts = int(machine.get_symmetry_part())
        self.content += ['-- parts      = {}'.format(parts)]
        if geom.is_stator():
            if machine.get_num_slots() > 0:
                self.content += [
                    '-- num_slots  = {}'.format(
                        machine.get_num_slots())
                ]
                # fix stator subregions
                sregs = [a.type for a in geom.list_of_areas()
                         if a.type in (TYPE_YOKE, TYPE_TOOTH)]
                if len(sregs) > 0 and set(sregs) == {TYPE_TOOTH}:
                    for a in geom.list_of_areas():
                        if a.type == TYPE_TOOTH:
                            a.type = TYPE_YOKE
                            logger.debug("FIXED STATOR SUBREGION")

        self.content += ['-- min_radius = {}'.format(geom.min_radius),
                         '-- max_radius = {}'.format(geom.max_radius),
                         '-- min_corner = {}, {}'.format(
                             geom.start_min_corner(0),
                             geom.start_min_corner(1)),
                         '-- max_corner start = {}, {}'.format(
                             geom.start_max_corner(0),
                             geom.start_max_corner(1)),
                         '\n']
        if inner:
            self.content += ['inner_da_start = {}'
                             .format(geom.dist_start_max_corner()),
                             'inner_da_end = {}'
                             .format(geom.dist_end_max_corner())]
        if outer:
            slice = 1 if machine.is_mirrored() else 2
            self.content += [
                '-- create air layer outside',
                'x0, y0 = {}, {}'.format(
                    geom.start_max_corner(0),
                    geom.start_max_corner(1)),
                'hair = 1.0',
                'parts = {}'.format(machine.get_num_parts()),
                'rmax = {}'.format(geom.max_radius),
                'r1 = rmax + hair',
                'r, phi = c2pr(x0, y0)',
                'x1, y1 = pr2c(r1, phi)',
                'x2, y2 = pr2c(r1, {}*math.pi/parts)'.format(slice),
                f'-- end max corner {geom.end_corners[-1]}',
                f'-- center {geom.center}',
                'r2 = {}'.format(geom.dist_end_max_corner(mirrored=False)),
                'x3, y3 = pr2c(r2, {}*math.pi/parts)'.format(slice),
                'nc_line(x0, y0, x1, y1, 0)',
                'nc_circle_m(x1, y1, x2, y2, 0.0, 0.0, 0)',
                'nc_line(x2, y2, x3, y3, 0)',
                'x0, y0 = pr2c(rmax + hair/2, math.pi/parts/2)',
                'create_mesh_se(x0, y0)',
                '\n',
                'outer_da_start = {}'.format(
                    geom.dist_start_min_corner()),
                'outer_da_end = {}'.format(
                    geom.dist_end_min_corner())
            ]
        if self.mtype == 'PMSM':
            self.content += [
                         'xmag = {}',
                         'ymag = {}',
                         'mag_orient = {}']
        self.content += ['\n',
                         'mag_exists = 0',
                         'if mcvkey_yoke == nil then',
                         '  mcvkey_yoke = "dummy"',
                         '  ur = 1000.0',
                         'end',
                         'x0_iron_tooth, y0_iron_tooth = 0.0, 0.0',
                         'x0_iron_yoke, y0_iron_yoke = 0.0, 0.0',
                         'x0_shaft, y0_shaft = 0.0, 0.0',
                         '\n']

        subregions = {}
        num_windings = 0
        num_magnets = 0
        magor = []

        for area in geom.list_of_areas():
            if area.number_of_elements() > 1:
                p = area.get_point_inside(geom)
                if p:
                    self.content.append("x0, y0 = {}, {}".format(p[0], p[1]))
                    # self.content.append("point(x0, y0, red, 4)")  # for debugging
                    self.content.append("create_mesh_se(x0, y0)")

                if area.is_winding():
                    if area.type not in subregions:
                        subregions[area.type] = 1
                    if self.mtype == 'PMSM' or outer:
                        num_windings += 1
                        rmin, rmax = area.minmax_dist_from_center((0,0))
                        self.content.append(
                            f'rcoil_{num_windings} = {rmin}, {rmax}')
                        self.content.append('m.xcoil_{}, m.ycoil_{} = x0, y0'.
                                            format(num_windings, num_windings))
                    else:
                        self.content.append('m.xcoil_r, m.ycoil_r = x0, y0')

                elif area.is_magnet():
                    if area.type not in subregions:
                        subregions[area.type] = 1
                    num_magnets += 1
                    self.content.append('xmag[{}], ymag[{}] = x0, y0'.
                                        format(num_magnets, num_magnets))
                    self.content.append('mag_orient[{}] = {}'.
                                        format(num_magnets, area.phi))
                    magor.append([p[0], p[1], area.phi])  # must correct rotation?

                elif area.type > 0:
                    if area.type in subregions:
                        self.content.append(
                            'add_to_subreg(x0, y0, "{}")'.
                            format(area.name()))
                    else:
                        subregions[area.type] = 1
                        color = area.color()
                        if color == 'cyan':
                            color = 'skyblue'
                        self.content.append(
                            'def_new_subreg(x0, y0, "{}", "{}")'.
                            format(area.name(), color))
                    if area.is_stator_iron_yoke():
                        self.content.append(
                            'x0_iron_yoke, y0_iron_yoke = x0, y0')
                    elif area.is_stator_iron_tooth():
                        self.content.append(
                            'x0_iron_tooth, y0_iron_tooth = x0, y0')
                    elif area.is_rotor_iron():
                        self.content.append(
                            'x0_iron_yoke, y0_iron_yoke = x0, y0')
                    elif area.is_shaft():
                        self.content.append(
                            'x0_shaft, y0_shaft = x0, y0')

                self.content.append("\n")

        txt = ["if x0_iron_yoke > 0.0 then",
               "  if mcvkey_yoke ~= 'dummy' then",
               '    def_mat_fm_nlin(x0_iron_yoke, y0_iron_yoke, "blue", mcvkey_yoke, 100)',
               '  else',
               '    def_mat_fm(x0_iron_yoke, y0_iron_yoke, ur, 100)',
               '  end',
               'end\n']
        self.content.append('\n'.join(txt))

        txt = ["if x0_iron_tooth > 0.0 then",
               "  if(x0_iron_yoke == 0 and mcvkey_yoke ~= 'dummy') then",
               "    def_mat_fm_nlin(x0_iron_tooth, y0_iron_tooth, 'blue', mcvkey_yoke, 100)",
               "  else",
               "    if (mcvkey_teeth ~= 'dummy' and mcvkey_teeth ~= nil) then",
               "      def_mat_fm_nlin(x0_iron_tooth, y0_iron_tooth, 'blue', mcvkey_teeth, 100)",
               "    else",
               "      def_mat_fm(x0_iron_tooth, y0_iron_tooth, ur, 100)",
               "    end",
               "  end",
               'end\n']
        self.content.append('\n'.join(txt))

        txt = ['if x0_shaft > 0.0 then',
               "  if mcvkey_shaft ~= 'dummy' then",
               '    def_mat_fm_nlin(x0_shaft, y0_shaft, "lightgrey", mcvkey_shaft, 100)',
               '  else',
               '    def_mat_fm(x0_shaft, y0_shaft, ur, 100)',
               '  end',
               'end\n']
        self.content.append('\n'.join(txt))

        if num_windings > 0:
            if geom.is_mirrored():
                self.content += [
                    'r, phi = c2pr(m.xcoil_1, m.ycoil_1)',
                    'm.xcoil_2, m.ycoil_2 = pr2c(r, {}*2.0 - phi)'.format(geom.alfa)]
            self.content.append('m.wdg_location  = 1 -- stator\n')

        if num_magnets > 0:
            self.content.append('mag_exists   = {}'.
                                format(num_magnets))
            if geom.is_mirrored():
                self.content.append('mag_mirrored = true')
            else:
                self.content.append('mag_mirrored = false')
            self.content.append('mag_alpha    = {}\n'
                                .format(geom.get_alfa()))

        if geom.is_mirrored():
            self.content.append('-- mirror')
            self.content.append('mirror_nodechains({}, {}, {}, {})\n'.format(
                geom.mirror_corners[1][0],   # max x1
                geom.mirror_corners[1][1],   # max y1
                geom.mirror_corners[0][0],   # min x2
                geom.mirror_corners[0][1]))  # min y2

        self.content.append('parts_gen = {}'.format(geom.num_variable()))

        # angle after mirroring
        self.content.append('alfa = {}\n'.format(geom.get_alfa()))

        self.content += ['-- rotate',
                         'x1, y1 = {}, {}'.format(
                             geom.start_corners[0][0],
                             geom.start_corners[0][1])]  # min xy1
        if outer:
            self.content.append('x2, y2 = pr2c(r1, 0.0)')
        else:
            self.content.append('x2, y2 = {}, {}'.format(
                geom.start_corners[1][0],
                geom.start_corners[1][1]))  # max xy1

        if geom.is_mirrored():
            self.content.append('x3, y3 = pr2c(x2, alfa)')
            self.content.append('x4, y4 = pr2c(x1, alfa)')
        elif outer:
            self.content += ['x3, y3 = pr2c(r1, 2*math.pi/parts+phi)',
                             'x4, y4 = {}, {}'.format(
                                 geom.end_corners[0][0],
                                 geom.end_corners[0][1])]  # min xy4
        else:
            self.content += ['x3, y3 = {}, {}'.format(
                                 geom.end_corners[-1][0],
                                 geom.end_corners[-1][1]), # min xy3
                             'x4, y4 = {}, {}'.format(
                                 geom.end_corners[0][0],
                                 geom.end_corners[0][1]),  # min xy4
                             'def_bcond_vpo(x4, y4, x1, y1)']

        self.content.append('if parts_gen > 1 then')
        if geom.corners_dont_match():
            txt = ['  -- Warning: corners dont match',
                   '  mirror_nodechains(x3, y3, x4, y4)',
                   '  if {} > 2 then'.format(geom.num_variable()),
                   '    my_alfa = alfa * 2',
                   '    x3, y3 = pr2c(x2, my_alfa)',
                   '    x4, y4 = pr2c(x1, my_alfa)',
                   '    n = {} / 2'.format(geom.num_variable()),
                   '    rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,n-1)',
                   '  end']
            self.content.append('\n'.join(txt))
        else:
            self.content.append(
                '  rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,parts_gen-1)')
        self.content.append('end')

        self.content.append('alfa = parts_gen * alfa\n')

        if self.fm_nlin:
            self.content.append('\nx0, y0 = {}, {}'. format(
                self.fm_nlin[0], self.fm_nlin[1]))

            mat = ["if fm_nlin_mcvfile ~= 'dummy' then",
                   "  if fm_nlin_mcvfile == 'air' then",
                   "    def_mat_fm(x0,y0, 1.0, fm_nlin_rlen)",
                   "  else",
                   "    def_mat_fm_nlin(x0,y0, fm_nlin_colour, fm_nlin_mcvfile, fm_nlin_rlen)",
                   "  end",
                   "else",
                   "  def_mat_fm(x0,y0, 1000.0, fm_nlin_rlen)",
                   "end"]
            self.content.append('\n'.join(mat))

        if self.shaft:
            mat = ['\nif shaft_mat==1 then',
                   '  def_mat_fm_nlin(0.1,0.1,"lightgrey",fm_nlin_mcvfile_shft,fm_nlin_rlen)',
                   'end']
            self.content.append('\n'.join(mat))

        if num_magnets:
            # TODO: this is clearly a hack to adapt the mag orientation
            # for cases where magnets have multiple small slices
            rotangle = 0
            mag = np.array(magor)
            x, y = mag.T[:2]
            if np.allclose(mag.T[2] - np.arctan2(y, x), np.pi/2, atol=2e-2):
                # all angles are 90°
                rotangle = 90
            self.content += [
                '-- pm magnets',
                'if mag_exists > 0 then',
                '  if( m.remanenc == nil ) then',
                '     m.remanenc = 1.15',
                '     m.relperm  = 1.05',
                '  end',
                '  if( m.orient == nil ) then',
                '    m.orient = m.parallel',
                '  end',
                '  if( m.magncond == nil ) then',
                '    m.magncond = 625000.0',
                '  end',
                '  if( m.rlen == nil ) then',
                '    m.rlen = 100',
                '  end',
                '  for i=0, m.npols_gen-1 do',
                '    for n=1, mag_exists do',
                '      r, p = c2pr(xmag[n], ymag[n])',
                '      x0, y0 = pr2c(r, i*mag_alpha+p)',
                '      phi = (i*mag_alpha+mag_orient[n])*180/math.pi',
                '      if ( i % 2 == 0 ) then',
                f'        phi = phi - {rotangle}',
                '        color = "red"',
                '      else',
                f'        phi = phi - {rotangle} + 180',
                '        color = "green"',
                '      end',
                '      if(m.mcvkey_magnet == nil) then',
                '        def_mat_pm(x0, y0, color, m.remanenc, m.relperm,',
                '                   phi, m.orient, m.magncond, m.rlen)',
                '      else',
                '        def_mat_pm_nlin(x0, y0, color, m.mcvkey_magnet,',
                '                 phi, m.orient, m.magncond, m.rlen)',
                '      end',
                '      if mag_mirrored then',
                '        x0, y0 = pr2c(r, (i+1)*mag_alpha-p)',
                '        phi = ((i+1)*mag_alpha-mag_orient[n])*180/math.pi',
                '        if ( i % 2 == 0 ) then',
                f'          phi = phi - {rotangle}',
                '        else',
                f'          phi = phi - {rotangle} + 180',
                '        end',
                '        if(m.mcvkey_magnet == nil) then',
                '          def_mat_pm(x0, y0, color, m.remanenc, m.relperm,',
                '                       phi, m.orient, m.magncond, m.rlen)',
                '        else',
                '          def_mat_pm_nlin(x0, y0, color, m.mcvkey_magnet,',
                '                   phi, m.orient, m.magncond, m.rlen)',
                '        end',
                '      end',
                '    end',
                '  end',
                'end']

        return self.content

    def render_main(self,
                    m_inner, m_outer,
                    inner, outer,
                    params):
        '''create fsl statements for main file'''
        if not (m_inner and m_outer):
            logger.warning("ERROR: Rotor or Stator missing")
            return []

        if params['external_rotor']:
            stator = inner
            rotor = outer
        else:
            stator = outer
            rotor = inner

        num_layers = min(m_inner.num_of_layers() +
                         m_outer.num_of_layers(),
                         2)
        self.agndst = params.get('agndst', 0.1)
        excwin = []
        if self.mtype == 'EESM':
            excwin = [
                'r,beta  = c2pr(m.xcoil_r, m.ycoil_r)',
                'phi = math.pi/m.num_poles',
                'alpha = phi-beta',
                'm.num_wires = 1',
                'dir = {"wi", "wo"}',
                'xcoil, ycoil = pr2c(r, phi - alpha)',
                'def_new_wdg(xcoil, ycoil, "violet", "Exc", m.num_wires, 10.0, dir[1])',
                'xcoil, ycoil = pr2c(r, phi + alpha)',
                'add_to_wdg(xcoil, ycoil, wsamekey, dir[2], "wser")',
                'for i = 2, m.npols_gen do',
                '  n = (i+1) % 2 + 1',
                '  phi = phi + 2*math.pi/m.num_poles',
                '  xcoil, ycoil = pr2c(r, phi - alpha)',
                '  add_to_wdg(xcoil, ycoil, wsamekey, dir[n], "wser")',
                '  xcoil, ycoil = pr2c(r, phi + alpha)',
                '  n = i % 2 + 1',
                '  add_to_wdg(xcoil, ycoil, wsamekey, dir[n], "wser")',
                'end']
        return [
            '-- generated from DXF by femagtools {}'.format(__version__),
            'exit_on_error = false',
            'exit_on_end = false',
            'verbosity = 2\n',
            'new_model_force("{}","Test")'.format(self.model),
            'global_unit(mm)',
            'pickdist(0.001)',
            'cosys(polar)\n',
            'dy1 = {}'.format(params.get('dy1', 0.0)),
            'da1 = {}'.format(params.get('da1', 0.0)),
            'dy2 = {}'.format(params.get('dy2', 0.0)),
            'da2 = {}'.format(params.get('da2', 0.0)),
            'ag  = (da1 - da2)/2\n',
            'm.tot_num_slot   = {}'.format(params.get('tot_num_slot', 0)),
            'm.num_sl_gen     = {}'.format(params.get('num_sl_gen', 0)),
            'm.num_poles      = {}'.format(params.get('num_poles', 0)),
            'm.num_pol_pair   = m.num_poles/2',
            'm.num_slots      = m.num_sl_gen',
            'm.npols_gen      = m.num_poles * m.num_sl_gen / m.tot_num_slot',
            'm.tot_num_sl     = m.tot_num_slot',
            'm.fc_radius      = (da1+da2)/4',
            'm.fc_radius1     = m.fc_radius',
            'm.arm_length     = 1.0',
            'pre_models("basic_modpar")\n',
            'm.airgap         = 2*ag/3',
            'm.nodedist       = 1.0',
            'agndst           = {}'.format(self.agndst),
            "mcvkey_teeth = 'dummy'",
            "mcvkey_yoke = 'dummy'",
            "mcvkey_shaft = 'dummy'",
            "ur = 1000.0",
            "ndt(agndst)"] + stator + [
                "mcvkey_teeth = 'dummy'",
                "mcvkey_yoke = 'dummy'",
                "mcvkey_shaft = 'dummy'",
                "ndt(agndst)"] + rotor + [
                    '-- airgap',
                    'ndt(agndst)',
                    'r1 = m.fc_radius - ag/6',
                    'x1, y1 = pr2c(r1, alfa)',
                    'n = math.floor(m.fc_radius*alfa/agndst + 1.5)',
                    'nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)\n',
                    'r2 = m.fc_radius + ag/6',
                    'x2, y2 = pr2c(r2, alfa)',
                    'nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)\n',
                    'if inner_da_start == nil then',
                    '  inner_da_start = da2/2',
                    'end',
                    'x1, y1 = pr2c(inner_da_start, 0.0)',
                    'nc_line(x1, y1, r1, 0.0, 0.0)\n',
                    'if outer_da_start == nil then',
                    '  outer_da_start = da1/2',
                    'end',
                    'x2, y2 = pr2c(outer_da_start, 0.0)',
                    'nc_line(r2, 0.0, x2, y2, 0.0)\n',
                    'if m.tot_num_slot > m.num_sl_gen then',
                    '  if inner_da_end == nil then',
                    '    inner_da_end = inner_da_start',
                    '  end',
                    '  x3, y3 = pr2c(inner_da_end, alfa)',
                    '  x4, y4 = pr2c(r1, alfa)',
                    '  nc_line(x3, y3, x4, y4, 0, 0)\n',
                    '  if outer_da_end == nil then',
                    '    outer_da_end = outer_da_start',
                    '  end',
                    '  x3, y3 = pr2c(outer_da_end, alfa)',
                    '  x4, y4 = pr2c(r2, alfa)',
                    '  nc_line(x3, y3, x4, y4, 0, 0)',
                    'end\n',
                    'x0, y0 = pr2c(r1-ag/6, alfa/2)',
                    'create_mesh_se(x0, y0)',
                    'x0, y0 = pr2c(r2+ag/6, alfa/2)',
                    'create_mesh_se(x0, y0)\n',
                    'connect_models()\n',
                    '-- Gen_winding',
                    'if m.xcoil_1 ~= nil then',
                    '  m.num_phases      = 3',
                    '  m.num_layers      = {}'.format(num_layers),
                    '  m.num_wires       = 1',
                    '  m.coil_span       = {}'.format(
                        params.get(
                            'tot_num_slot', 1)//params.get('num_poles', 1)),
                    '  m.current         = 0.0',
                    '  m.mat_type        = 1.0 -- rotating',
                    '  m.wind_type       = 1.0 -- winding & current',
                    '  m.win_asym        = 1.0 -- sym',
                    '  m.wdg_location    = 1.0 -- stator',
                    '  m.curr_inp        = 0.0 -- const',
                    '  m.dq_offset       = 0',
                    '  pre_models("Gen_winding")',
                    '  pre_models("gen_pocfile")',
                    'end\n',
                    'save_model(cont)\n'] + excwin
