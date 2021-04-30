# -*- coding: utf-8 -*-
"""
    femagtools.dxfsl.plotrenderer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    renders geometry to plot

"""
import numpy as np
from .geom import convex_hull
import logging
try:
    import matplotlib
    import matplotlib.pyplot as pl
    import matplotlib.patches as pch

    matplotlibversion = (int(matplotlib.__version__.split('.')[0]) +
                         int(matplotlib.__version__.split('.')[1])/10)
except ImportError:  # ModuleNotFoundError:
    matplotlibversion = 0  # no matplotlib

logger = logging.getLogger(__name__)


class PlotRenderer(object):
    """renders geom with matplotlib"""
    def __init__(self):
        self.fig = None
        self.background = '#eeeeee'
#        self.background = facecolor='lightblue',
        if matplotlibversion == 0:
            import matplotlib   # throws exception
        pass

    def figure(self):
        if self.fig is None:
            self.fig = pl.figure(figsize=(9, 10))
            #         facecolor='lightblue'
            if matplotlibversion > 2:
                pl.tight_layout(h_pad=0.2, w_pad=0.2)
            pl.subplots_adjust(bottom=0.05, top=0.95,
                               hspace=0.25, wspace=0.15)
        return self.fig

    def add_plot(self, rows=1, cols=1, num=1):
        self.ax = self.figure().add_subplot(rows, cols,
                                            num)
#                                            facecolor=self.background)
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

    def ellipse(self, center, width, height,
                rtheta, start_param, end_param, color='blue'):
        self.ax.add_patch(pch.Arc(center,
                                  width,
                                  height,
                                  angle=rtheta*180/np.pi,
                                  theta1=start_param*180/np.pi,
                                  theta2=end_param*180/np.pi,
                                  color=color))

    def line(self, p1, p2, color='blue', e=None):
        self.ax.add_line(pl.Line2D((p1[0], p2[0]),
                                   (p1[1], p2[1]), color=color))

    def circle(self, center, radius, color='blue'):
        self.ax.add_patch(pch.Circle(center, radius,
                                     fill=False, color=color))

    def point(self, p, marker, color='blue'):
        pl.plot([p[0]], [p[1]], marker, color=color)

    def text(self, p, txt):
        pl.text(p[0], p[1], txt)

    def fill(self, x, y, color, alpha):
        self.ax.fill(x, y, color, alpha=alpha, edgecolor='blue')

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
                        c = s.center_of_connection(4)
                        pl.plot([c[0]], [c[1]], 'gs')

        geom.remove_areas(incl_bnd)

        # draw all remaining edges and circles
        [circle.render(self) for circle in geom.circles()]
        [attr['object'].render(self)
         for e1, e2, attr in geom.g.edges(data=True)]
        if draw_center:
            for e1, e2, attr in geom.g.edges(data=True):
                c = attr['object'].center_of_connection(4)
                pl.plot([c[0]], [c[1]], 'gs')

        geom.render_cut_lines(self)
        self.ax.axis('scaled')
        self.ax.set_aspect('equal')
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
            self.point(p, 'o', 'green')
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
        write_id = kwargs.get('write_id', False)
        title = kwargs.get('title', "")
        show = kwargs.get('show', True)
        write_png = kwargs.get('png', False)
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
                    self.ax.axis('scaled')
                    self.ax.set_aspect('equal')
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
                                            num)
#                                            facecolor=self.background)
        if len(title) > 0:
            self.ax.set_title(title, size=14)
        self.ax.grid(color='blue', linewidth=0.5)

        if with_hull:
            for h in convex_hull(geom.virtual_nodes()):
                pl.plot([h[0]], [h[1]], 'ro')

        if with_corners:
            for c in geom.start_corners:
                self.point(c, 's')
            for c in geom.end_corners:
                self.point(c, 's')

        if draw_center:
            for area in geom.list_of_areas():
                for s in area.elements():
                    c = s.center_of_connection(4)
                    pl.plot([c[0]], [c[1]], 'gs')

        if draw_inside:
            for area in geom.list_of_areas():
                p = area.get_point_inside(geom)
                if p:
                    pl.plot([p[0]], [p[1]], 'o', color='magenta')
                    if write_id:
                        self.text(p, area.identifier())

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
                self.point(p, 'o', color='red')

        if with_center and geom.center:
            self.point(geom.center, 'o', color='darkgreen')
            x_min = min(x_min, geom.center[0]-5)
            x_max = max(x_max, geom.center[0]+5)
            y_min = min(y_min, geom.center[1]-5)
            y_max = max(y_max, geom.center[1]+5)

        self.ax.axis('scaled')
        self.ax.set_aspect('equal')
        
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

        if show:
            if write_png:
                self.write_plot(title)
            else:
                self.show_plot()

    def render_areas(self, geom, **kwargs):
        with_nodes = kwargs.get('with_nodes', False)
        single_view = kwargs.get('single_view', False)
        title = kwargs.get('title', "geometry")
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
        a = len(geom.list_of_areas())
        n = 0
        for area in geom.list_of_areas():
            n += 1
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
                    self.point(p, 'o', color='magenta')

                if single_view:
                    self.ax.axis('scaled')
                    self.ax.set_aspect('equal')

                    mytitle = "Area {} of {} in {}".format(n, a, title)
                    self.ax.set_title(mytitle, size=14)
                    pl.show()

        if not single_view:
            self.ax.axis('scaled')
            self.ax.set_aspect('equal')
            if show:
                self.show_plot()

    def render_area(self, area):
        fig = pl.figure()
        self.ax = fig.add_subplot(111)

        area.render(self, 'red', with_nodes=True)

        self.ax.axis('scaled')
        self.ax.set_aspect('equal')       
        pl.show()

    def draw_slot(self, id, slot, ax):
        poly = pch.Polygon(slot, fill=True, color='#ee82ee')

        center = (sum(list(zip(*slot))[0])/len(slot),
                  sum(list(zip(*slot))[1])/len(slot))
        ax.text(center[0], center[1], str(id))
        ax.add_patch(poly)

    def write_plot(self, name):
        filename = '{}.png'.format(name)
        fig = self.figure()
        fig.savefig(filename)
