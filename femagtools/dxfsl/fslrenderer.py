# -*- coding: utf-8 -*-
"""
    dxfsl
    ~~~~~

    renders geometry


"""
import numpy as np
from .shape import Shape
import logging
import io

logger = logging.getLogger(__name__)


#############################
#        FslRenderer        #
#############################

class FslRenderer(object):
    """a model that can created by FSL"""
    header = ['exit_on_error = false',
              'exit_on_end = false',
              'verbosity = 2',
              'pickdist = {}',
              'fm_nlin_colour=blue',
              'fm_nlin_color=fm_nlin_colour',
              'fm_nlin_mcvfile="M470P65A-AM-50Hz-96"',
              'fm_nlin_mcvfile_shft="ST52-3"',
              'fm_nlin_rlen = 100',
              'shaft_mat=1',
              'agndst = 0.5',
              'new_model_force("{}","Test")',
              'blow_up_wind(0,0, {}, {})']

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

    def render(self, machine, filename, inner=False, outer=False):
        '''create file with nodechains'''

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
        self.content.append(u'-- parts      = {}'
                            .format(parts))
        self.content.append(u'-- min_radius = {}'
                            .format(geom.min_radius))
        self.content.append(u'-- max_radius = {}'
                            .format(geom.max_radius))
        self.content.append(u'-- min_corner = {}, {}'
                            .format(geom.start_min_corner(0),
                                    geom.start_min_corner(1)))
        self.content.append(u'-- max_corner = {}, {}'
                            .format(geom.start_max_corner(0),
                                    geom.start_max_corner(1)))
        self.content.append(u'\n')

        if inner:
            self.content.append(u'tmp.inner_max_corner_x = {}'
                                .format(geom.start_max_corner(0)))
            self.content.append(u'tmp.inner_max_corner_y = {}'
                                .format(geom.start_max_corner(1)))
        if outer:
            self.content.append(u'tmp.outer_min_corner_x = {}'
                                .format(geom.start_min_corner(0)))
            self.content.append(u'tmp.outer_min_corner_y = {}'
                                .format(geom.start_min_corner(1)))
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
                    self.content.append(u'tmp.xmag[{}], tmp.ymag[{}] = x0, y0'.
                                        format(num_magnets, num_magnets))
                    self.content.append(u'tmp.mag_orient[{}] = {}'.
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
            self.content.append(u'tmp.coil_exists   = {}'.
                                format(num_windings))
            if geom.is_mirrored():
                self.content.append(u'tmp.coil_mirrored = true')
            else:
                self.content.append(u'tmp.coil_mirrored = false')
            self.content.append(u'tmp.coil_alpha    = {}'.format(geom.alfa))
            self.content.append(u'm.wdg_location  = 1.0 -- stator\n')

        if num_magnets > 0:
            self.content.append(u'tmp.mag_exists   = {}'.
                                format(num_magnets))
            if geom.is_mirrored():
                self.content.append(u'tmp.mag_mirrored = true')
            else:
                self.content.append(u'tmp.mag_mirrored = false')
            self.content.append(u'tmp.mag_alpha    = {}\n'
                                .format(geom.get_alfa()))

        if geom.is_mirrored():
            self.content.append(u'-- mirror')
            self.content.append(u'mirror_nodechains({}, {}, {}, {})\n'.format(
                geom.mirror_corners[1][0],   # max x1
                geom.mirror_corners[1][1],   # max y1
                geom.mirror_corners[0][0],   # min x2
                geom.mirror_corners[0][1]))  # min y2

        # Winkel nach allfälligem Spiegeln
        self.content.append(u'alfa = {}\n'.format(geom.get_alfa()))

        self.content.append(u'-- rotate')
        self.content.append(u'x1, y1 = {}, {}'.format(
                geom.start_corners[0][0], geom.start_corners[0][1]))  # min xy1
        self.content.append(u'x2, y2 = {}, {}'.format(
                geom.start_corners[1][0], geom.start_corners[1][1]))  # max xy2

        if geom.is_mirrored():
            self.content.append(u'x3, y3 = pr2c(x2, alfa)')
            self.content.append(u'x4, y4 = pr2c(x1, alfa)')
        else:
            self.content.append(u'x3, y3 = {}, {}'.format(
                    geom.end_corners[1][0], geom.end_corners[1][1]))  # max xy3
            self.content.append(u'x4, y4 = {}, {}'.format(
                    geom.end_corners[0][0], geom.end_corners[0][1]))  # min xy4

        self.content.append(u'if {} > 1 then'.format(geom.num_variable()))
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

        # f.write(u"\nadapt_window()\n")
        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write('\n'.join(self.content))

    def render_main(self,
                    motor,
                    m_inner,
                    m_outer,
                    params,
                    filename):
        '''create main file'''
        if not (m_inner and m_outer):
            logger.warning("ERROR: Rotor or Stator missing")
            return

        self.content = []
        self.content.append(u'exit_on_error = false')
        self.content.append(u'exit_on_end = false')
        self.content.append(u'verbosity = 2\n')

        self.content.append(u'new_model_force("{}","Test")'.
                            format(self.model))
        self.content.append(u'global_unit(mm)')
        self.content.append(u'pickdist(0.001)')
        self.content.append(u'cosys(polar)\n')

        geom_inner = m_inner.geom
        geom_outer = m_outer.geom

        self.content.append(u'tmp = {}')
        self.content.append(u'tmp.xmag = {}')
        self.content.append(u'tmp.ymag = {}')
        self.content.append(u'tmp.mag_orient = {}')
        self.content.append(u'tmp.mag_exists = 0')
        self.content.append(u'tmp.coil_exists = 0\n')

        self.content.append(u'dy1 = {}'.format(params.get('dy1', 0.0)))
        self.content.append(u'da1 = {}'.format(params.get('da1', 0.0)))
        self.content.append(u'dy2 = {}'.format(params.get('dy2', 0.0)))
        self.content.append(u'da2 = {}'.format(params.get('da2', 0.0)))
        self.content.append(u'ag  = (da1 - da2)/2\n')

        self.content.append(u'm.tot_num_slot   = {}'
                            .format(params.get('tot_num_slot', 0)))
        self.content.append(u'm.num_sl_gen     = {}'
                            .format(params.get('num_sl_gen', 0)))
        self.content.append(u'm.num_poles      = {}'
                            .format(params.get('num_poles', 0)))
        self.content.append(u'm.num_pol_pair   = m.num_poles/2')
        self.content.append(u'm.num_slots      = m.num_sl_gen')
        self.content.append(u'm.npols_gen      = m.num_poles * m.num_sl_gen / m.tot_num_slot')
        self.content.append(u'm.tot_num_sl     = m.tot_num_slot')
        self.content.append(u'm.fc_radius      = (da1+da2)/4')
        self.content.append(u'm.fc_radius1     = m.fc_radius')
        self.content.append(u'm.arm_length     = 1.0')
        self.content.append(u'pre_models("basic_modpar")\n')

        self.content.append(u'm.airgap         = 2*ag/3')
        self.content.append(u'm.nodedist       = 1.0')
        # set of useful node angles
        dagset = [np.pi/45, np.pi/90, np.pi/180, np.pi/360, np.pi/720]
        r = (params.get('da1', 0.0) + params.get('da2', 0.0))/4
        ag = abs((params.get('da1', 0.0) - params.get('da2', 0.0)))/6
        i = np.argmin(np.abs(np.array(dagset) - np.arctan2(ag, r)))
        self.content.append(u'agndst           = {}'.format(dagset[i-1]*r))

        self.content.append(u'blow_up_wind(0, 0, 10, 10)\n')

        self.content.append(u"mcvkey_yoke = 'dummy'")
        self.content.append(u"mcvkey_shaft = 'dummy'")
        self.content.append(u"ur = 1000.0")
        self.content.append(u"ndt(agndst)")
        self.content.append(u'dofile("{}_{}.fsl")\n'
                            .format(self.model, geom_inner.kind))

        self.content.append(u"mcvkey_yoke = 'dummy'")
        self.content.append(u"mcvkey_shaft = 'dummy'")
        self.content.append(u"ndt(agndst)")
        self.content.append(u'dofile("{}_{}.fsl")\n'
                            .format(self.model, geom_outer.kind))

        # Airgap
        txt = [u'-- airgap',
               u'ndt(agndst)',
               u'r1 = da2/2 + ag/3',
               u'x1, y1 = pr2c(r1, alfa)',
               u'n = r1*alfa/agndst + 1',
               u'nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)\n',
               u'r2 = da2/2 + 2*ag/3',
               u'x2, y2 = pr2c(r2, alfa)',
               u'nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)\n',
               u'if tmp.inner_max_corner_x == nil then',
               u'  tmp.inner_max_corner_x = da2/2',
               u'end',
               u'x1, y1 = tmp.inner_max_corner_x, 0.0',
               u'nc_line(x1, y1, r1, 0.0, 0.0)\n',
               u'if tmp.outer_min_corner_x == nil then',
               u'  tmp.outer_min_corner_x = da1/2',
               u'end',
               u'x2, y2 = tmp.outer_min_corner_x, 0.0',
               u'nc_line(r2, 0.0, x2, y2, 0.0)\n',
               u'x3, y3 = pr2c(x1, alfa)',
               u'x4, y4 = pr2c(r1, alfa)',
               u'nc_line(x3, y3, x4, y4, 0, 0)\n',
               u'x3, y3 = pr2c(x2, alfa)',
               u'x4, y4 = pr2c(r2, alfa)',
               u'nc_line(x3, y3, x4, y4, 0, 0)\n',
               u'x0, y0 = pr2c(r1-ag/6, alfa/2)',
               u'create_mesh_se(x0, y0)',
               u'x0, y0 = pr2c(r2+ag/6, alfa/2)',
               u'create_mesh_se(x0, y0)\n']
        self.content.append(u'\n'.join(txt))

        self.content.append(u'connect_models()\n')

        # Windings
        txt = [u'-- Gen_winding',
               u'if tmp.coil_exists > 0 then',
               u'  m.num_phases      = 3']
        self.content.append(u'\n'.join(txt))

        num_layers = min(geom_inner.num_of_windings() +
                         geom_outer.num_of_windings(),
                         2)
        self.content.append(u'  m.num_layers      = {}'
                            .format(num_layers))
        txt = [u'  m.num_wires       = 1',
               u'  m.coil_span       = 1',
               u'  m.current         = 0.0',
               u'  m.mat_type        = 1.0 -- rotating',
               u'  m.wind_type       = 1.0 -- winding & current',
               u'  m.win_asym        = 1.0 -- sym',
               u'  m.wdg_location    = 1.0 -- stator',
               u'  m.curr_inp        = 0.0 -- const',
               u'  m.dq_offset       = 0',
               u'  if tmp.coil_exists == 1 and tmp.coil_mirrored then',
               u'    r, phi = c2pr(m.xcoil_1, m.ycoil_1)',
               u'    m.xcoil_2, m.ycoil_2 = pr2c(r, tmp.coil_alpha*2.0 - phi)',
               u'  end\n',
               u'  pre_models("Gen_winding")',
               u'end\n']
        self.content.append(u'\n'.join(txt))

        txt = [u'-- pm magnets',
               u'if tmp.mag_exists > 0 then',
               u'  alfa = tmp.mag_alpha',
               u'  m.remanenc = 1.15',
               u'  m.relperm  = 1.05',
               u'  for i=0, m.npols_gen-1 do',
               u'    for n=1, tmp.mag_exists do',
               u'      r, p = c2pr(tmp.xmag[n], tmp.ymag[n])',
               u'      phi = i*alfa+p',
               u'      x0, y0 = pr2c(r, phi)',
               u'      phi_orient = i*alfa+tmp.mag_orient[n]',
               u'      orient = phi_orient*180/math.pi',
               u'      if ( i % 2 == 0 ) then',
               u'        def_mat_pm(x0, y0, red, m.remanenc, m.relperm,',
               u'                   orient-180, m.parallel, 100)',
               u'      else',
               u'        def_mat_pm(x0, y0, green, m.remanenc, m.relperm,',
               u'                   orient, m.parallel, 100)',
               u'      end',
               u'      if tmp.mag_mirrored then',
               u'        phi = (i+1)*alfa-p',
               u'        x0, y0 = pr2c(r, phi)',
               u'        phi_orient = (i+1)*alfa-tmp.mag_orient[n]',
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
               u'end',
               u'save_model(cont)']
        self.content.append(u'\n'.join(txt))

        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write('\n'.join(self.content))
