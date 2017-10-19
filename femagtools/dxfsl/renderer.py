# -*- coding: utf-8 -*-
"""
    dxfsl
    ~~~~~

    renders geometry

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.

"""
import numpy as np
import numpy.linalg as la

import matplotlib.pylab as pl
import matplotlib.patches as pch
import femagtools.dxfsl.geom as g
import logging


logger = logging.getLogger(__name__)


def get_point_inside(area):
    """return point inside area"""
    center = [n.center_of_connection() for n in area]
    xymin = np.min(center, axis=0)
    xymax = np.max(center, axis=0)

    dx, dy = xymax - xymin
    logger.debug("path min %s max %s dx, %f dy %f", xymin, xymax, dx, dy)

    # intersect path with a line somewhere in the center
    if dx > dy:
        line = g.Line(g.Element(start=(xymin[0] + dx/2, xymin[1]),
                                end=(xymin[0] + dx/2, xymax[1])))
    else:
        line = g.Line(g.Element(start=(xymin[0], xymin[1] + dy/2),
                                end=(xymax[0], xymin[1] + dy/2)))

    # collect all intersecting edges with line
    intersect = []
    pickdist = 0.1
    logger.debug("line %s -- %s", line.start(), line.end())
    for e in area:
        intersect += [(round(ip[0], 2), round(ip[1], 2))
                      for ip in e.intersect(line, pickdist,
                                            include_end=True) if ip]

    if len(set(intersect)) > 1:
        logger.debug("Intersections %s", list(set(intersect)))
        # take the first 2 intersection points
        p0, p1 = np.array(sorted(set(intersect))[:2])
        return (p0 + (p1-p0)/2.).tolist()
    return ()


class PlotRenderer(object):
    def __init__(self):
        pass
    
    def arc(self, startangle, endangle, center, radius, color = 'blue'):
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
        
    def render(self, geom, filename=None, **kwargs):
        draw_center = kwargs.get('draw_center', False)
        draw_corners = kwargs.get('draw_corners', False)
        draw_hull = kwargs.get('draw_hull', False)
        incl_bnd = kwargs.get('incl_bnd', False)

        fig = pl.figure()
        self.ax = fig.add_subplot(111)

        if draw_hull:
            for h in g.convex_hull(geom.g.nodes()):
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

#                p = get_point_inside(area)
#                if p:
#                    pl.plot([p[0]], [p[1]], 'ro')
#                    self.ax.text(p[0], p[1], str(id))

        geom.remove_areas(incl_bnd)

        # draw all remaining edges and circles
        [circle.render(self) for circle in geom.circles()]
        [attr['object'].render(self)
         for e1, e2, attr in geom.g.edges(data=True)]
        if draw_center:
            for e1, e2, attr in geom.g.edges(data=True):
                c = attr['object'].center_of_connection()
                pl.plot([c[0]], [c[1]], 'gs')

        if draw_corners:
            for c in geom.find_corners(geom.g):
                pl.plot([c[0]], [c[1]], 'bs')

        geom.render_schnitt(self)
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
                x = lower.pop()
#                self.point(x, 'ro')
            lower.append(p)
        print("lower={}".format(lower))

    def render_elements(self, geom, type, **kwargs):
        with_nodes = kwargs.get('with_nodes', False)
        with_hull = kwargs.get('with_hull', False)
        with_corners = kwargs.get('with_corners', False)

        fig = pl.figure()
        self.ax = fig.add_subplot(111)
        
        for e in geom.elements(type):
            e.render(self, 'blue', with_nodes)
        
        if with_hull:
            for h in g.convex_hull(geom.virtual_nodes()):
                pl.plot([h[0]], [h[1]], 'ro')

        if with_corners:
#            for c in geom.find_corners(geom.g.nodes()):
#                self.point(c, 'rs')
            for c in g.find_corners(geom.g.nodes(), True):
                self.point(c, 'rs')

        geom.render_schnitt(self)
        geom.render_airgaps(self)

        if geom.center:
            self.circle(geom.center, 3, 'darkgreen')
            
        self.ax.axis('scaled', aspect='equal')
        pl.show()

    def render_areas(self, geom, **kwargs):
        with_nodes = kwargs.get('with_nodes', False)
        
        fig = pl.figure()
        self.ax = fig.add_subplot(111)

        colors = ('red', 'green', 'blue', 'magenta', 'orange', 'grey', 'darkgreen')
        
        c = -1
        for area in geom.areas(True):
            if len(area) > 1:
                c += 1
                if c >= len(colors):
                    c = 0
                for s in area:
                    s.render(self, colors[c], with_nodes)

        self.ax.axis('scaled', aspect='equal')
        pl.show()
            
    def draw_slot(self, id, slot, ax):
        poly = pch.Polygon(slot, fill=True, color='#ee82ee')

        center = (sum(list(zip(*slot))[0])/len(slot),
                  sum(list(zip(*slot))[1])/len(slot))
        ax.text(center[0], center[1], str(id))
        ax.add_patch(poly)

        
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
                if not e in handled:
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

    def circle(self, center, radius):
        num = int(2*np.pi*radius)
        if num < 8:
            num = 8
        circle = [u'cx, cy = {}, {}'.format(center[0],
                                            center[1]),
                  u'nc_circle_m(cx,cy,{}, {})'.format(radius, num),
                  u'create_mesh_se(cx, cy)\n']
        self.content += circle

    def arc(self, startangle, endangle, center, radius):
        num = 0
        if self.nodedist > 0:
            s = startangle
            d = endangle - s
            n = int(d/(2*np.pi))
            if n < 0:
                n -= 1
            d -= n*2*np.pi
            num = int(radius*d/self.nodedist + 1)
            if num < 3 and radius*d > self.nodedist:
                num = 3

        p1 = (center[0] + radius*np.cos(startangle),
              center[1] + radius*np.sin(startangle))
        p2 = (center[0] + radius*np.cos(endangle),
              center[1] + radius*np.sin(endangle))
        self.content.append(
            u"nc_circle_m({}, {}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1], num))

    def line(self, p1, p2):
        num = 0
        if self.nodedist > 0:
            l = la.norm(np.asarray(p1)-p2)
            num = int(l/self.nodedist + 1)
        self.content.append(
            u"nc_line({}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1], num))

    def set_nodedist(self, e, airgapOutside):
        """determine position of e and adjust nodedist if needed"""
        r = la.norm(e.center_of_connection())
        n = -1
        for level in self.ndttab:
            if ((airgapOutside and level[0] < r) or
                (not airgapOutside and level[0] > r)):
                break
            n += 1
        if n != self.ndtpos:
            self.ndtpos = n
            self.content.append(
                u"ndt({})\n".format(
                    self.ndttab[n][1] if n >= 0 else "agndst"))

    def render(self, geom, airgapOutside=True, incl_bnd=True):
        '''create fsl commands with nodechains'''

        handled = set()  # prevent multiple creation
        if airgapOutside:
            dy, da = geom.diameters[0], geom.diameters[-1]
            self.ndttab = [(0.98*da/2, '1.6*agndst'),
                           (dy/2 + 0.4*(da-dy), '3*agndst')]
        else:
            dy, da = geom.diameters[-1], geom.diameters[0]
            self.ndttab = [(1.06*da/2, 1.6), (dy/2 - 0.2*(da-da), 3)]
        self.content = []
        self.ndtpos = -1
        self.content.append(u'\n\nndt(agndst)\n')
            
        # collect all areas
        coll = []
        for i, area in enumerate(geom.areas(incl_bnd)):
            if len(area) > 1:
                r = la.norm(np.sum([p.center_of_connection()
                                    for p in area], axis=0)/len(area))
                coll.append((r, area))
        geom.remove_areas(incl_bnd)
        # and all remaining edges and circles
        for circle in geom.circles():
            coll.append((la.norm(circle.center), circle))

        for e1, e2, attr in geom.g.edges(data=True):
            r = la.norm(attr['object'].center_of_connection())
            coll.append((r, attr['object']))

        # create nodechains in sorted order and set nodedist
        self.nodedist = 0
        for r, e in sorted(coll, key=lambda x: x[0], reverse=airgapOutside):
            if isinstance(e, list):
                self.content.append(u"-- {})\n".format(r))
                acoll = sorted([(la.norm(p.center_of_connection()), p)
                                for p in e],
                               key=lambda x: x[0],
                               reverse=airgapOutside)
                for p in acoll:
                    if p in handled:
                        continue
                    self.set_nodedist(p[1], airgapOutside)
                    p[1].render(self)
                    handled.add(p)
                    
                p = get_point_inside(e)
                self.content.append(
                    u"create_mesh_se({}, {})\n".format(p[0], p[1]))

            else:
                self.set_nodedist(e, airgapOutside)
                e.render(self)

        if airgapOutside:
            self.content.append(u'\nx0, y0 = pr2c({}, {})'. format(
                1.01*dy/2, geom.alpha/2))
        else:
            self.content.append(u'\nx0, y0 = pr2c({}, {})'. format(
                0.99*dy/2, geom.alpha/2))
        if not incl_bnd:
            self.content.append(u'create_mesh_se(x0, y0)')
        self.content.append(u'def_new_subreg(x0, y0, "Rotor", green)')

        if geom.mirror_axis:
            self.content.append(u'\nmirror_nodechains({})\n'.format(
                ', '.join([str(x) for x in geom.mirror_axis])))

        self.content += [u'-- rotate',
                         u'alfa = {}'.format(2*geom.alpha),
                         u'x1, y1 = {}, 0.0'.format(
                             dy/2 if airgapOutside else da/2),
                         u'x2, y2 = {}, 0.0'.format(
                             da/2 if airgapOutside else dy/2),
                         u'x3, y3 = pr2c(x2, alfa)',
                         u'x4, y4 = pr2c(x1, alfa)',
                         u'rotate_copy_nodechains(x1,y1,x2,y2,x3,y3,x4,y4,ncopies)']
        
        mat = [u"if mcvkey_yoke ~= 'dummy' then",
               u"  if fm_nlin_mcvfile == 'air' then",
               u"    def_mat_fm(x0,y0, 1.0, fm_nlin_rlen)",
               u"  else",
               u"    def_mat_fm_nlin(x0,y0, blue, mcvkey_yoke, fm_nlin_rlen)",
               u"  end",
               u"else",
               u"  def_mat_fm(x0,y0, 1000.0, fm_nlin_rlen)",
               u"end"]
        self.content += mat

        return self.content
