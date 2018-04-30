# -*- coding: utf-8 -*-
"""
    dxfsl
    ~~~~~

    renders geometry


"""
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
import matplotlib.patches as pch
from .shape import Shape
from .geom import convex_hull
import logging
import io

logger = logging.getLogger(__name__)
matplotlibversion = int(matplotlib.__version__.split('.')[0])
matplotlibrelease = int(matplotlib.__version__.split('.')[1])


#############################
#       PlotRenderer        #
#############################

class PlotRenderer(object):
    def __init__(self):
        self.fig = None
        self.background = '#eeeeee'
        pass

    def figure(self):
        if self.fig is None:
            self.fig = pl.figure(facecolor='lightblue',
                                 figsize=(9, 10))
            if matplotlibversion == 2 and matplotlibrelease > 0:
                pl.tight_layout(h_pad=0.2, w_pad=0.2)
            pl.subplots_adjust(bottom=0.05, top=0.95, hspace=0.25, wspace=0.15)
        return self.fig

    def add_plot(self, rows=1, cols=1, num=1):
        self.ax = self.figure().add_subplot(rows, cols,
                                            num,
                                            facecolor=self.background)
        return self.ax

    def add_emptyplot(self, rows=1, cols=1, num=1, title=''):
        ax = self.add_plot(rows, cols, num)
        ax.axis('off')
        if title:
            ax.set_title(title, size=14)

    def new_legend_handle(self, color, alpha, text):
        return pch.Patch(color=color, alpha=alpha, label=text)

    def show_plot(self):
        if self.fig is not None:
            pl.show()
            self.fig = None

    def arc(self, startangle, endangle, center, radius, color='blue'):
        self.ax.add_patch(pch.Arc(center, 2*radius, 2*radius,
                                  angle=0,
                                  theta1=startangle*180/np.pi,
                                  theta2=endangle*180/np.pi,
                                  color=color))

    def line(self, p1, p2, color='blue'):
        self.ax.add_line(pl.Line2D((p1[0], p2[0]),
                                   (p1[1], p2[1]), color=color))

    def circle(self, center, radius, color='blue'):
        self.ax.add_patch(pch.Circle(center, radius,
                                     fill=False, color=color))

    def point(self, p, col, color='blue'):
        pl.plot([p[0]], [p[1]], col, color=color)

    def fill(self, x, y, color, alpha):
        self.ax.fill(x, y, 'b', alpha=alpha, color=color, edgecolor='blue')

    def fill_circle(self, center, radius, color, alpha):
        circle = pl.Circle(center, radius, color=color, alpha=alpha)
        self.ax.add_artist(circle)
        circle = pl.Circle(center, radius, color='blue', fill=False)
        self.ax.add_artist(circle)

    def render(self, geom, filename=None, **kwargs):
        draw_center = kwargs.get('draw_center', False)
        draw_hull = kwargs.get('draw_hull', False)
        incl_bnd = kwargs.get('incl_bnd', False)

        fig = pl.figure()
        self.ax = fig.add_subplot(111)

        if draw_hull:
            for h in convex_hull(geom.g.nodes()):
                pl.plot([h[0]], [h[1]], 'ro')

        id = 0
        for area in geom.areas(incl_bnd):
            id += 1
            logger.info("area %d (size %d)", id, len(area))
            if len(area) > 1:
                for s in area:
                    s.render(self)
                    if draw_center:
                        c = s.center_of_connection()
                        pl.plot([c[0]], [c[1]], 'gs')

        geom.remove_areas(incl_bnd)

        # draw all remaining edges and circles
        [circle.render(self) for circle in geom.circles()]
        [attr['object'].render(self)
         for e1, e2, attr in geom.g.edges(data=True)]
        if draw_center:
            for e1, e2, attr in geom.g.edges(data=True):
                c = attr['object'].center_of_connection()
                pl.plot([c[0]], [c[1]], 'gs')

        geom.render_cut_lines(self)
        self.ax.axis('scaled', aspect='equal')
        if filename:
            pl.savefig(filename)
        else:
            pl.show()

    def print_convex_hull(self, nodes):
        points = list(nodes)
        points.sort()
        if len(points) <= 1:
            return points

        def cross(o, a, b):
            dx = a[0] - o[0], a[1] - o[1]
            dy = b[0] - o[0], b[1] - o[1]
            return dx[0] * dy[1] - dx[1] * dy[0]

        print("=== Print Convex Hull ===")
        lower = []
        for p in points:
            self.point(p, 'go')
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)
        print("lower={}".format(lower))

    def render_elements(self, geom, type, **kwargs):
        with_nodes = kwargs.get('with_nodes', False)
        with_hull = kwargs.get('with_hull', False)
        with_corners = kwargs.get('with_corners', False)
        with_center = kwargs.get('with_center', False)
        single_view = kwargs.get('single_view', False)
        neighbors = kwargs.get('neighbors', False)
        draw_center = kwargs.get('draw_center', False)
        draw_inside = kwargs.get('draw_inside', False)
        fill_areas = kwargs.get('fill_areas', False)
        title = kwargs.get('title', "")
        show = kwargs.get('show', True)
        rows = kwargs.get('rows', 1)
        cols = kwargs.get('cols', 1)
        num = kwargs.get('num', 1)
        points = kwargs.get('points', [])

        if show:
            rows = 1
            cols = 1
            num = 1

        mm = geom.minmax()
        x_min, x_max = mm[0]-5, mm[1]+5
        y_min, y_max = mm[2]-5, mm[3]+5

        if single_view:
            count = 0
            for e in geom.elements(type):
                print("Render Element {}".format(e))
                if count == 0:
                    self.ax = self.figure().add_subplot(111)
                    self.ax.axis('scaled', aspect='equal')
                    self.ax.set_xlim(x_min, x_max)
                    self.ax.set_ylim(y_min, y_max)
                e.render(self, 'blue', True)

                count += 1
                if count == 3:
                    pl.show()
                    count = 0

            if count != 0:
                self.show_plot()
            return

        self.ax = self.figure().add_subplot(rows, cols,
                                            num,
                                            facecolor=self.background)
        if len(title) > 0:
            self.ax.set_title(title, size=14)
        self.ax.grid(color='blue', linewidth=0.5)

        if with_hull:
            for h in convex_hull(geom.virtual_nodes()):
                pl.plot([h[0]], [h[1]], 'ro')

        if with_corners:
            for c in geom.start_corners:
                self.point(c, 'rs')
            for c in geom.end_corners:
                self.point(c, 'rs')

        if draw_center:
            for area in geom.list_of_areas():
                for s in area.elements():
                    c = s.center_of_connection()
                    pl.plot([c[0]], [c[1]], 'gs')

        if draw_inside:
            for area in geom.list_of_areas():
                p = area.get_point_inside(geom)
                if p:
                    pl.plot([p[0]], [p[1]], 'ro', color='magenta')

        if fill_areas:
            handles = geom.render_area_fill(self)
            if handles:
                legend = pl.legend(handles=handles, loc='best',
                                   fancybox=True,
                                   framealpha=1.0,
                                   fontsize=8,
                                   labelspacing=0.2)
                frame = legend.get_frame()
                frame.set_facecolor('white')
                frame.set_edgecolor('blue')

        for e in geom.elements(type):
            e.render(self, 'blue', with_nodes)

        geom.render_cut_lines(self)
        geom.render_airgaps(self)
        if neighbors:
            geom.render_neighbors(self)

        if len(points) > 0:
            for p in points:
                self.point(p, 'ro', color='red')

        if with_center and geom.center:
            self.point(geom.center, 'ro', color='darkgreen')
            x_min = min(x_min, geom.center[0]-5)
            x_max = max(x_max, geom.center[0]+5)
            y_min = min(y_min, geom.center[1]-5)
            y_max = max(y_max, geom.center[1]+5)

        self.ax.axis('scaled', aspect='equal')
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

        if show:
            self.show_plot()

    def render_areas(self, geom, **kwargs):
        with_nodes = kwargs.get('with_nodes', False)
        single_view = kwargs.get('single_view', False)
        title = kwargs.get('title', "")
        show = kwargs.get('show', True)
        rows = kwargs.get('rows', 1)
        cols = kwargs.get('cols', 1)
        num = kwargs.get('num', 1)

        if show:
            rows = 1
            cols = 1
            num = 1

        if not single_view:
            self.ax = self.figure().add_subplot(rows, cols, num)
            if len(title) > 0:
                self.ax.set_title(title, size=14)

        colors = ('red', 'green', 'blue', 'magenta',
                  'orange', 'grey', 'darkgreen')

        c = -1
        for area in geom.list_of_areas():
            if area.number_of_elements() > 1:
                c += 1
                if c >= len(colors):
                    c = 0
                if single_view:
                    print("*** AREA has {} elements ***".
                          format(area.number_of_elements()))
                    fig = pl.figure()
                    self.ax = fig.add_subplot(111)

                area.render(self, colors[c], with_nodes)
                p = area.get_point_inside(geom)
                if p:
                    self.point(p, 'ro', color='magenta')

                if single_view:
                    self.ax.axis('scaled', aspect='equal')
                    pl.show()

        if not single_view:
            self.ax.axis('scaled', aspect='equal')
            if show:
                self.show_plot()

    def render_area(self, area):
        fig = pl.figure()
        self.ax = fig.add_subplot(111)

        area.render(self, 'red', with_nodes=True)

        self.ax.axis('scaled', aspect='equal')
        pl.show()

    def draw_slot(self, id, slot, ax):
        poly = pch.Polygon(slot, fill=True, color='#ee82ee')

        center = (sum(list(zip(*slot))[0])/len(slot),
                  sum(list(zip(*slot))[1])/len(slot))
        ax.text(center[0], center[1], str(id))
        ax.add_patch(poly)


#############################
#       DumpRenderer        #
#############################

class DumpRenderer(object):
    def __init__(self, name):
        self.model = name

    def circle(self, center, radius):
        self.content.append(u'Crcle({}, {}, {})'.format(center[0],
                                                        center[1],
                                                        radius))

    def arc(self, startangle, endangle, center, radius):
        p1 = (center[0] + radius*np.cos(startangle),
              center[1] + radius*np.sin(startangle))
        p2 = (center[0] + radius*np.cos(endangle),
              center[1] + radius*np.sin(endangle))
        self.content.append(
            u"Arc({}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1]))

    def line(self, p1, p2):
        self.content.append(
            u"Line({}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1]))

    def render(self, geom):
        '''create file with elements'''
        incl_bnd = True
        self.content = []
        handled = set()

        # collect all areas
        for i, area in enumerate(geom.areas(incl_bnd)):
            self.content.append('-- {}'.format(i))
            for e in area:
                if e not in handled:
                    e.render(self)
                    handled.add(e)
        self.content.append('--')
        geom.remove_areas(incl_bnd)
        # and all remaining edges and circles
        for circle in geom.circles():
            circle.render(self)

        for e1, e2, attr in geom.g.edges(data=True):
            attr['object'].render(self)

        return self.content


#############################
#       NewFslRenderer      #
#############################

class NewFslRenderer(object):
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

    def line(self, p1, p2, color='blue'):
        num = 0
#        if self.nodedist > 0:
#            l = la.norm(np.asarray(p1)-p2)
#            num = int(l/self.nodedist + 1)
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

    def render(self, geom, filename, inner=False, outer=False):
        '''create file with nodechains'''

        self.content = []

        self.content.append(u'\n\nndt(agndst)\n')

        ndt_list = [(0.25, 1.25), (0.5, 2), (0.75, 3.0), (1.1, 3.0)]
        dist = geom.max_radius - geom.min_radius
        el_sorted = self.sorted_elements(geom, inner)

        x = 0
        for d, e in el_sorted:
            d_percent = d / dist
            if ndt_list[x][0] < d_percent:
                self.content.append(u'\nndt({}*agndst)\n'.
                                    format(ndt_list[x][1]))
                while x < 3 and ndt_list[x][0] < d_percent:
                    x += 1
#            self.content.append(u'-- d={} / dist={} == {}'.
#                                format(d, dist, d_percent))
            e.render(self)

        self.content.append(u'\n')

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
                    self.content.append(u'm.xmag[{}], m.ymag[{}] = x0, y0'.
                                        format(num_magnets, num_magnets))
                    self.content.append(u'm.mag_orient[{}] = {}'.
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
                self.content.append(u"\n")

        if num_windings > 0:
            self.content.append(u'm.coil_exists   = {}'.
                                format(num_windings))
            if geom.is_mirrored():
                self.content.append(u'm.coil_mirrored = true')
            else:
                self.content.append(u'm.coil_mirrored = false')
            self.content.append(u'm.coil_alpha    = {}\n'.format(geom.alfa))

        if num_magnets > 0:
            self.content.append(u'm.mag_exists   = {}'.
                                format(num_magnets))
            if geom.is_mirrored():
                self.content.append(u'm.mag_mirrored = true')
            else:
                self.content.append(u'm.mag_mirrored = false')
            self.content.append(u'm.mag_alpha    = {}\n'.format(geom.get_alfa()))

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

        self.content.append(u'if m.{}_ncopies > 1 then'.format(geom.kind))
        self.content.append(
            u'  rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,m.{}_ncopies-1)'
            .format(geom.kind))
        self.content.append(u'end')

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

    def render_main(self, motor, m_inner, m_outer,
                    filename, with_header=False):
        '''create main file'''

        self.content = []
        self.content.append(u'exit_on_error = false')
        self.content.append(u'exit_on_end = false')
        self.content.append(u'verbosity = 2\n')

        self.content.append(u'new_model_force("{}","Test")'.
                            format(self.model))
        self.content.append(u'global_unit(mm)')
        self.content.append(u'pickdist(0.001)')
        self.content.append(u'cosys(polar)\n')

        geom_inner = None
        geom_outer = None

        if m_inner and m_outer:
            geom_inner = m_inner.geom
            geom_outer = m_outer.geom

            parts_inner = int(m_inner.get_symmetry_part())
            parts_outer = int(m_outer.get_symmetry_part())

            if parts_inner > parts_outer:
                num_slots = parts_inner
                num_poles = parts_outer
                npols_gen = int(geom_outer.get_symmetry_copies()+1)
            else:
                num_slots = parts_outer
                num_poles = parts_inner
                npols_gen = int(geom_inner.get_symmetry_copies()+1)
            self.content.append(u'm.xmag = {}')
            self.content.append(u'm.ymag = {}')
            self.content.append(u'm.mag_orient = {}')
            self.content.append(u'm.mag_exists = 0')
            self.content.append(u'm.coil_exists = 0\n')
            self.content.append(u'm.num_poles = {}'.format(num_poles))
            self.content.append(u'm.tot_num_slot = {}'.format(num_slots))

            self.content.append(u'm.npols_gen = {}'.format(npols_gen))
            self.content.append(
                u'm.num_slots = m.tot_num_slot*m.npols_gen/m.num_poles')

        self.content.append(u'da1 = {}'.format(
            2*geom_outer.min_radius))
        self.content.append(u'da2 = {}'.format(
            2*geom_inner.max_radius))
        self.content.append(u'ag = (da1 - da2)/2\n')

        if m_inner and m_outer:
            self.content.append(u'm.tot_num_sl  = m.tot_num_slot')
            self.content.append(u'm.fc_radius   = (da1+da2)/4')
            self.content.append(u'm.fc_radius1  = m.fc_radius')
            self.content.append(u'pre_models("basic_modpar")\n')

        self.content.append(
            u'agndst = math.pi*m.fc_radius*m.npols_gen/m.num_poles/90')

        if geom_inner:
            if parts_inner > parts_outer:
                ncopies = 'm.num_slots'
            else:
                ncopies = 'm.npols_gen'
            self.content.append(
                u'm.{}_ncopies = {}'.format(
                    geom_inner.kind, ncopies))

        if geom_outer:
            if parts_inner > parts_outer:
                ncopies = 'm.npols_gen'
            else:
                ncopies = 'm.num_slots'
            self.content.append(
                u'm.{}_ncopies = {}'.format(
                    geom_outer.kind, ncopies))

        self.content.append(u'blow_up_wind(0, 0, 10, 10)')

        if geom_inner:
            self.content.append(
                u'dofile("{}_{}.fsl")'.format(
                    self.model, geom_inner.kind))
        if geom_outer:
            self.content.append(
                u'dofile("{}_{}.fsl")'.format(
                    self.model, geom_outer.kind))

        alfa = geom_inner.get_alfa() * (geom_inner.get_symmetry_copies()+1)
        alfa2 = geom_outer.get_alfa() * (geom_outer.get_symmetry_copies()+1)
        assert(np.isclose(alfa, alfa2))

        self.content.append(u'alfa = {}\n'.format(alfa))

        # Airgap
        txt = [u'-- airgap',
               u'ndt(agndst)',
               u'r1 = da2/2 + ag/3',
               u'x1, y1 = pr2c(r1, alfa)',
               u'n = r1*alfa/agndst + 1',
               u'nc_circle_m(r1, 0, x1, y1, 0.0, 0.0, n)\n',
               u'r2 = da2/2 + 2*ag/3',
               u'x2, y2 = pr2c(r2, alfa)',
               u'nc_circle_m(r2, 0, x2, y2, 0.0, 0.0, n)\n']
        self.content.append(u'\n'.join(txt))
        self.content.append(u'x1, y1 = {}, {}'
                            .format(geom_inner.start_max_corner(0),
                                    geom_inner.start_max_corner(1)))
        self.content.append(u'nc_line(x1, y1, r1, 0.0, 0.0)\n')
        self.content.append(u'x2, y2 = {}, {}'
                            .format(geom_outer.start_min_corner(0),
                                    geom_outer.start_min_corner(1)))
        self.content.append(u'nc_line(r2, 0.0, x2, y2, 0.0)\n')

        txt = [u'x3, y3 = pr2c(x1, alfa)',
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
               u'if m.coil_exists > 0 then',
               u'  m.num_phases      = 3',
               u'  m.num_layers      = 2',
               u'  m.num_wires       = 1',
               u'  m.coil_span       = 1',
               u'  m.current         = 0.0',
               u'  m.mat_type        = 1.0 -- rotating',
               u'  m.wind_type       = 1.0 -- winding & current',
               u'  m.win_asym        = 1.0 -- sym',
               u'  m.wdg_location    = 1.0 -- stator',
               u'  m.curr_inp        = 0.0 -- const',
               u'  m.dq_offset       = 0',
               u'  if m.coil_exists == 1 and m.coil_mirrored then',
               u'    r, phi = c2pr(m.xcoil_1, m.ycoil_1)',
               u'    m.xcoil_2, m.ycoil_2 = pr2c(r, m.coil_alpha*2.0 - phi)',
               u'  end\n',
               u'  pre_models("Gen_winding")',
               u'end\n']
        self.content.append(u'\n'.join(txt))

        txt = [u'-- iron',
               u'urr    = 1000',
               u"mcvkey = 'dummy'"]
        self.content.append(u'\n'.join(txt))
        if geom_inner:
            points = geom_inner.get_points_in_iron()
            if points:
                self.content.append(u'x0 = {} -- {}'
                                    .format(points[0][0], geom_inner.kind))
                self.content.append(u'y0 = {}'.format(points[0][1]))
                self.content.append(u"if mcvkey ~= 'dummy' then")
                self.content.append(u'  def_mat_fm_nlin(x0, y0, blue, mcvkey, 100)')
                self.content.append(u'else')
                self.content.append(u'  def_mat_fm(x0, y0, 1000.0, 100)')
                self.content.append(u'end')
        if geom_outer:
            points = geom_outer.get_points_in_iron()
            if points:
                self.content.append(u'x0 = {} -- {}'
                                    .format(points[0][0], geom_outer.kind))
                self.content.append(u'y0 = {}'.format(points[0][1]))
                self.content.append(u"if mcvkey ~= 'dummy' then")
                self.content.append(u'  def_mat_fm_nlin(x0, y0, blue, mcvkey, 100)')
                self.content.append(u'else')
                self.content.append(u'  def_mat_fm(x0, y0, 1000.0, 100)')
                self.content.append(u'end')
        self.content.append(u'')

        txt = [u'-- pm magnets',
               u'if m.mag_exists > 0 then',
               u'  alfa = m.mag_alpha',
               u'  m.remanenc = 1.15',
               u'  m.relperm  = 1.05',
               u'  for i=0, m.npols_gen-1 do',
               u'    for n=1, m.mag_exists do',
               u'      r, p = c2pr(m.xmag[n], m.ymag[n])',
               u'      phi = i*alfa+p',
               u'      x0, y0 = pr2c(r, phi)',
               u'      phi_orient = i*alfa+m.mag_orient[n]',
               u'      orient = phi_orient*180/math.pi',
               u'      if ( i % 2 == 0 ) then',
               u'        def_mat_pm(x0, y0, red, m.remanenc, m.relperm,',
               u'                   orient-180, m.parallel, 100)',
               u'      else',
               u'        def_mat_pm(x0, y0, green, m.remanenc, m.relperm,',
               u'                   orient, m.parallel, 100)',
               u'      end',
               u'      if m.mag_mirrored then',
               u'        phi = (i+1)*alfa-p',
               u'        x0, y0 = pr2c(r, phi)',
               u'        phi_orient = (i+1)*alfa-m.mag_orient[n]',
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
        self.content.append(u'\n'.join(txt))

        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write('\n'.join(self.content))
