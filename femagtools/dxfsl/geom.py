# -*- coding: utf-8 -*-
"""
    NOTE: This code is in highly experimental state.
          Use at your own risk.

  Author: Ronald Tanner
    Date: 2017/07/06
"""
from __future__ import print_function

import dxfgrabber
import numpy as np
import networkx as nx
import logging
import sys
from .corner import Corner
from .area import Area
from .shape import Element, Shape, Circle, Arc, Line, Point
from .machine import Machine
from .functions import less_equal, less, greater
from .functions import distance, alpha_line, alpha_points, alpha_angle
from .functions import point, points_are_close, is_point_inside_region
from .functions import line_m, line_n
from .functions import middle_point_of_line, middle_point_of_arc
from .functions import middle_angle
from .functions import normalise_angle, is_same_angle
from .functions import part_of_circle, gcd
from .functions import point_in_region
import io

logger = logging.getLogger('femagtools.geom')


nxversion = int(nx.__version__.split('.')[0])
geom_mindist = 0.0


def plot_area(area):
    """plot area for debug purposes"""
    import matplotlib.pylab as pl
    import matplotlib.patches as pch
    fig = pl.figure()
    ax = fig.add_subplot(111)

    for e in area:
        if isinstance(e, Line):
            ax.add_line(pl.Line2D((e.p1[0], e.p2[0]),
                                  (e.p1[1], e.p2[1]),
                                  color='red'))
        elif isinstance(e, Arc):
            ax.add_patch(pch.Arc(e.center, 2*e.radius, 2*e.radius,
                                 angle=0,
                                 theta1=e.startangle*180/np.pi,
                                 theta2=e.endangle*180/np.pi,
                                 color='darkblue'))
        elif isinstance(e, Circle):
            ax.add_patch(pch.Circle(e.center,
                                    e.radius,
                                    fill=False, color='darkblue'))

    ax.axis('scaled', aspect='equal')
    pl.show()

#############################
#            geom           #
#############################

ndec = 6  # number of decimals to round to


def intersect_and_split(inp_elements, rtol, atol):
    logger.info("Load input elements ... ")
    out_elements = []
    for e in inp_elements:
        out_size = len(out_elements)
        intersect_and_split_element(e, out_elements, 0, out_size, rtol, atol)
    logger.info(" ... loaded")
    return out_elements


def intersect_and_split_element(el, out_elements, out_start,
                                out_size, rtol, atol):
    # appends splitted elements
    # Unchanged out_size prevents repeated processing in recursive calls
    for x in range(out_start, out_size):
        split_el = add_or_split(el, x, out_elements, rtol, atol)
        if len(split_el) > 0:
            for e in split_el:
                intersect_and_split_element(e, out_elements, x+1,
                                            out_size, rtol, atol)
            return
    out_elements.append(el)


def add_or_split(el, x, out_elements, rtol, atol):
    if out_elements[x] is None:
        return []
    points = el.intersect_shape(out_elements[x], rtol, atol, True)
    if len(points) > 0:
        split_elements = out_elements[x].split(points, rtol, atol)
        if len(split_elements) > 0:
            out_elements += split_elements
            out_elements[x] = None
        split_el = el.split(points, rtol, atol)
        return split_el

    return []


def add_or_join(geom, n1, n2, entity, rtol, atol):
    """ adds a new entity to graph or joins entity with existing
    geom: Geometry
    n1, n2: nodes
    entity
    """
    if n1 == n2:
        logger.debug("Tiny element with same node on both sides ignored: %s", n1)
    else:
        geom.add_edge(n1, n2, entity)


def get_nodes_of_paths(g, c):
    nodes = set()
    i = 1
    for s in c[:-1]:
        for e in c[i:]:
            logger.info("%s --%s:", s, e)
            for p in nx.all_simple_paths(g, s, e):
                logger.info("   %s", len(list(p)))
                nodes |= set([n for n in p if g.degree(n) < 3])
        i += 1
    return nodes


def convex_hull(nodes):
    """collect edges that build the convex hull."""

    # Get a local list copy of the nodes and sort them lexically.
    points = list(nodes)
    points.sort()
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of
    # their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        dx = a[0] - o[0], a[1] - o[1]
        dy = b[0] - o[0], b[1] - o[1]
        return dx[0] * dy[1] - dx[1] * dy[0]

    # Build lower hull
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the
    # beginning of the other list.
    return lower[:-1] + upper[:-1]


def ccw(a, b, c):
    """
    checks if 3 points are placed counter clockwise or in line
    """
    cp = np.cross(np.asarray(b)-np.asarray(a),
                  np.asarray(c)-np.asarray(a))
    if np.allclose(cp, (0.0,)):
        return 0
    if cp > 0:
        return 1
    return -1


def intersect(a, b):
    """checks if 2 lines intersect"""
    return ccw(a.p1, a.p2, b.p1) != ccw(a.p1, a.p2, b.p2) and \
        ccw(b.p1, b.p2, a.p1) != ccw(b.p1, b.p2, a.p2)


def polylines(entity):
    """returns a collection of bulged vertices
    http://www.afralisp.net/archive/lisp/Bulges1.htm
    """
    if isinstance(entity.points, list):
        points = [(p[0], p[1]) for p in entity.points]
    else:
        points = [(p[0], p[1]) for p in entity.points()]
    i = 0
    if isinstance(entity.vertices, list):
        vertices = entity.vertices
    else:
        vertices = entity.vertices()

    for v in vertices:
        if hasattr(v, 'bulge'):
            b = v.bulge
        else:
            b = v.get_dxf_attrib('bulge', 0.0)
        p1 = points[i]
        try:
            p2 = points[i+1]
        except:
            p2 = points[0]
        if b != 0.0:
            dx, dy = p2[0] - p1[0], p2[1] - p1[1]
            c = np.sqrt(dx**2 + dy**2)
            s = b * c/2
            r = ((c/2)**2 + s**2)/2/s
            g = np.arctan2(dx, dy)
            a = 2*np.arctan(b)
            pc = (p1[0] + r*np.sin(g - np.pi/2 + a),
                  p1[1] + r*np.cos(g - np.pi/2 + a))
            if b < 0:
                sa = np.arctan2(p2[1]-pc[1], p2[0]-pc[0])
                ea = np.arctan2(p1[1]-pc[1], p1[0]-pc[0])
            else:
                sa = np.arctan2(p1[1]-pc[1], p1[0]-pc[0])
                ea = np.arctan2(p2[1]-pc[1], p2[0]-pc[0])
            logger.debug("Poly p1 %s p2 %s r %s", p1, p2, r)
            yield Arc(Element(center=pc, radius=np.abs(r),
                              start_angle=sa/np.pi*180,
                              end_angle=ea/np.pi*180))
        else:
            yield Line(Element(start=p1, end=p2))
        i += 1


def lw_polyline(entity):
    """returns a collection of bulged vertices
    http://www.afralisp.net/archive/lisp/Bulges1.htm
    """
    if isinstance(entity.points, list):
        points = [(p[0], p[1]) for p in entity.points]
    else:
        points = [(p[0], p[1]) for p in entity.points()]

    if points:
        p1 = points[0]
        for p2 in points[1:]:
            yield Line(Element(start=p1, end=p2))
            p1 = p2
    if entity.is_closed:
        yield Line(Element(start=p1, end=points[0]))


def spline(entity, min_dist=0.001):
    if False:
        yield Line(Element(start=entity.control_points[0],
                           end=entity.control_points[-1]))
        return

    if False:
        p_prev = None
        for p in entity.control_points:
            if p_prev:
                yield Line(Element(start=p_prev, end=p))
            p_prev = p
        return

    points_between = entity.control_points[1:-1]
    p1 = entity.control_points[0]
    pe = entity.control_points[-1]
    for p2 in points_between:
        dist_12 = distance(p1, p2)
        dist_2e = distance(p2, pe)
        if dist_2e < min_dist:
            logger.debug("SPLINE: ignor small end-distance {}".format(dist_2e))
            yield Line(Element(start=p1, end=pe))
            return

        if dist_12 > min_dist:
            yield Line(Element(start=p1, end=p2))
            p1 = p2
        else:
            logger.debug("SPLINE: ignor small distance {}".format(dist_12))

    yield Line(Element(start=p1, end=pe))


def insert_block(insert_entity, block, min_dist=0.001):
    logger.info('Insert {} entities from block {}'
                .format(len(block), insert_entity.name))

    logger.debug('Insert = {}'.format(insert_entity.insert))
    logger.debug('Rotation = {}'.format(insert_entity.rotation))
    logger.debug('Scale = {}'.format(insert_entity.scale))
    logger.debug('Rows = {}'.format(insert_entity.row_count))
    logger.debug('Cols = {}'.format(insert_entity.col_count))
    logger.debug('Row spacing = {}'.format(insert_entity.row_spacing))
    logger.debug('Col spacing = {}'.format(insert_entity.col_spacing))

    if insert_entity.insert != (0.0, 0.0, 0.0):
        logger.error('Different Location in Insert not supported')
        return

    if not (insert_entity.scale == (1.0, 1.0, 1.0) or
            insert_entity.scale == (1.0, 1.0, 0.0)):
        logger.error('Block scaling in Insert not supported')
        logger.error('  scale = {}'.format(insert_entity.scale))
        return

    if(insert_entity.row_count > 1 or
       insert_entity.col_count > 1 or
       insert_entity.row_spacing > 0 or
       insert_entity.col_spacing > 0):
        logger.error('Multi Block references in Insert not supported')
        return

    if insert_entity.rotation != 0.0:
        logger.error('Block Insert with rotation not supported')
        return

    for e in block:
        if e.dxftype == 'ARC':
            yield Arc(e)
        elif e.dxftype == 'CIRCLE':
            logger.debug("Circle %s, Radius %f", e.center[:2], e.radius)
            yield Circle(e)
        elif e.dxftype == 'LINE':
            yield Line(e)
        elif e.dxftype == 'POLYLINE':
            for p in polylines(e):
                yield p
        elif e.dxftype == 'LWPOLYLINE':
            for p in lw_polyline(e):
                yield p
        elif e.dxftype == 'SPLINE':
            for l in spline(e, min_dist=min_dist):
                yield l
        elif e.dxftype == 'INSERT':
            logger.warn("Nested Insert of Blocks not supported")
        else:
            logger.warn("unknown type %s in block %s",
                        e.dxftype, insert_entity.name)


def dxfshapes0(dxffile, mindist=0.01, layers=[]):
    """returns a collection of dxf entities (ezdxf)"""
    import ezdxf
    dwg = ezdxf.readfile(dxffile)
    id = 0
    # $ACADVER: AC1006 = R10, AC1009 = R11 and R12, AC1012 = R13,
    #   AC1014 = R14 AC1015 = Release 2000/0i/2
    # check units:
    # dwg.header['$ANGDIR'] 1 = Clockwise angles, 0 = Counterclockwise
    # dwg.header['$AUNITS'] Decimal Degrees, Deg/Min/Sec, Grads, Radians
    # dwg.header['$INSUNIT'] 1 = Inches; 2 = Feet; 3 = Miles;
    #   4 = Millimeters; 5 = Centimeters; 6 = Meters
    # dwg.header['$LUNITS']
    for e in dwg.modelspace():
        if e.dxftype() == 'ARC':
            yield Arc(e.dxf)
        elif e.dxftype() == 'CIRCLE':
            logger.debug("Circle %s, Radius %f", e.center[:2], e.radius)
            yield Circle(e.dxf)
        elif e.dxftype() == 'LINE':
            yield Line(e.dxf)
        elif e.dxftype() == 'POLYLINE':
            for p in polylines(e):
                yield p
        elif e.dxftype() == 'SPLINE':
            for l in spline(e, min_dist=mindist):
                yield l
        else:
            logger.info("Id %d4: unknown type %s", id, e.dxftype)
        id += 1


def dxfshapes(dxffile, mindist=0.01, layers=[]):
    """returns a collection of dxf entities (dxfgrabber)"""
    dwg = dxfgrabber.readfile(dxffile)
    # print("Layers = {}".format(dwg.layers.names()))
    id = 0
    # $ACADVER: AC1006 = R10, AC1009 = R11 and R12, AC1012 = R13,
    #   AC1014 = R14 AC1015 = Release 2000/0i/2
    # check units:
    # dwg.header['$ANGDIR'] 1 = Clockwise angles, 0 = Counterclockwise
    # dwg.header['$AUNITS'] Decimal Degrees, Deg/Min/Sec, Grads, Radians
    # dwg.header['$INSUNIT'] 1 = Inches; 2 = Feet; 3 = Miles;
    #   4 = Millimeters; 5 = Centimeters; 6 = Meters
    # dwg.header['$LUNITS']
    for e in dwg.modelspace():
        if not layers or e.layer in layers:
            if e.dxftype == 'ARC':
                yield Arc(e)
            elif e.dxftype == 'CIRCLE':
                logger.debug("Circle %s, Radius %f", e.center[:2], e.radius)
                yield Circle(e)
            elif e.dxftype == 'LINE':
                yield Line(e)
            elif e.dxftype == 'POLYLINE':
                for p in polylines(e):
                    yield p
            elif e.dxftype == 'LWPOLYLINE':
                for p in lw_polyline(e):
                    yield p
            elif e.dxftype == 'SPLINE':
                for l in spline(e, min_dist=mindist):
                    yield l
            elif e.dxftype == 'INSERT':
                block = dwg.blocks[e.name]
                for l in insert_block(e, block, min_dist=mindist):
                    yield l
            else:
                logger.info("Id %d4: unknown type %s", id, e.dxftype)
            id += 1


def is_connected(a, b):
    if a == b:
        return False
    return a[0] == b[0] or a[1] == b[0] or a[0] == b[1] or a[1] == b[1]


def single_path(edges):
    """sort edges to create a sequence of connected edges"""
    if edges:
        p = [edges[0]]
        remaining = edges[1:]
        while remaining:
            edges = remaining[:]
            remaining = []
            for e in edges:
                if is_connected(p[-1], e):
                    p.append(e)
                elif is_connected(p[0], e):
                    p.insert(0, e)
                else:
                    remaining.append(e)
        return p
    return []


#############################
#         Geometrie         #
#############################

class Geometry(object):
    """collection of connected shapes"""
    def __init__(self, elements=[],
                 rtol=1e-03,
                 atol=1e-03,
                 split=False,
                 debug=False):
        self._name = ''
        self.kind = ''
        self.mirror_corners = []
        self.start_corners = []
        self.end_corners = []
        self.sym_part = 0
        self.sym_counterpart = 0
        self.alfa = 0.0
        self.center = []
        self.min_radius = 0.0
        self.max_radius = 0.0
        self.is_inner = False
        self.is_outer = False
        self.cut_lines = []
        self.sym_area = None
        self.airgaps = []
        self.area_list = []
        self.g = nx.Graph()
        self.rtol = rtol
        self.atol = atol
        self.debug = debug
        self.num_edges = 0
        i = 0

        def get_elements(elements, split):
            if split:
                return intersect_and_split(elements, self.rtol, self.atol)
            else:
                return elements

        for e in get_elements(elements, split):
            if e:
                e.id = i
                n = self.find_nodes(e.start(), e.end())
                try:
                    add_or_join(self, n[0], n[1], e, self.rtol, self.atol)
                except Exception as ex:
                    logger.warn("EXCEPTION %s", ex)
                    if e:  # must be a circle
                        self.g.add_node(e.center, object=e)
            i += 1

        self.num_edges = self.number_of_edges()

    def shaft(self):
        """returns shaft diameter if any"""
        radius = []
        for a in self.arcs():
            if np.isclose(a.center, 0.0).all():
                logger.info("Shaft candidate (%f, %f) %r",
                            a.startangle/np.pi*180,
                            a.endangle/np.pi*180,
                            a.radius)
                if 0 < a.radius < self.diameters[-1]/2:
                    radius.append(a.radius)
        return 2*sorted(radius)[0] if radius else 0

    def set_diameters(self):
        try:
            min, max = self.min(), self.max()
            logger.info("Min %s max %s", min, max)
            self.diameters = (round(2*np.sqrt(min[0]**2 + min[1]**2), 3),
                              round(2*np.sqrt(max[0]**2 + max[1]**2), 3))
        except:
            self.diameters = []

    def move(self, offset):
        """moves all objects"""
        for e in self.g.edges(data=True):
            e[2]['object'].move(offset)
        for c in self.circles():
            c.move((-offset[0], -offset[1]))
        mapping = {n: (round(n[0]+offset[0], 3), round(n[1]+offset[1], 3))
                   for n in self.g.nodes()}
        nx.relabel_nodes(self.g, mapping, copy=False)

    def rotate(self, alpha):
        """rotates all objects by angle alpha"""
        T = np.array(((np.cos(alpha), -np.sin(alpha)),
                      (np.sin(alpha), np.cos(alpha))))
        for e in self.g.edges(data=True):
            e[2]['object'].transform(T)
        for c in self.circles():
            c.transform(T)
        rotnodes = np.dot(T, np.asarray(self.g.nodes()).T).T.tolist()
        mapping = {n: (round(r[0], 3), round(r[1], 3))
                   for n, r in zip(self.g.nodes(), rotnodes)}
        nx.relabel_nodes(self.g, mapping, copy=False)

    def scale(self, factor):
        """scales all objects"""
        for e in self.g.edges(data=True):
            e[2]['object'].scale(factor)
        for c in self.circles():
            c.scale(factor)
        mapping = {n: (round(factor*n[0], 3), round(factor*n[1], 3))
                   for n in self.g.nodes()}
        nx.relabel_nodes(self.g, mapping, copy=False)
        self.diameters = tuple([factor*d for d in self.diameters])

    def find_nodes0(self, *points):
        """return closest nodes to points in arg within pickdist"""
        return [tuple([round(int(x/self.atol+0.5)*self.atol, ndec)
                       for x in p])
                for p in points]

    def find_nodes(self, *points, **kwargs):
        """return closest nodes to points in arg within pickdist"""
        n = []
        nodes = list(kwargs.get('g', self.g))
        if nodes:
            anodes = np.asarray(nodes)
            for p in points:
                # la.norm on numpy below 1.8 does not accept axis
                c = anodes - p
                dist = np.sqrt(np.einsum('ij, ij->i', c, c))
                # dist = la.norm(np.asarray(nodes) - p, axis=1)
                idx = dist.argmin()
                if dist[idx] < self.atol:
                    n.append(nodes[idx])
                else:
                    n.append((round(p[0], ndec), round(p[1], ndec)))
        else:
            return [(round(p[0], ndec), round(p[1], ndec)) for p in points]
        return n

    def polyline_paths(self):
        """return lists of line paths"""
        paths = []
        g = self.g.copy()
        g.remove_edges_from(self.arcs())
        while g.number_of_edges() > 0:
            p = single_path(g.edges())
            g.remove_edges_from(p)
            # rearrange sequence to make it contiguous:
            px = [p[0]]
            for e in p[1:]:
                if px[-1][1] == e[1]:
                    px.append(tuple(reversed(e[:2])))
                elif px[0][0] == e[0]:
                    px.insert(0, tuple(reversed(e[:2])))
                elif px[0][0] == e[1]:
                    px.insert(0, e[:2])
                else:
                    px.append(e[:2])
            # only keep end point of each edge
            px[1:] = [x[-1] for x in px]
            px[0] = px[0][0]
            paths.append(px)
        return paths

    def number_of_nodes(self):
        """returns the number of nodes in graph"""
        return self.g.number_of_nodes()

    def number_of_edges(self):
        """return the number of edges in graph"""
        return self.g.number_of_edges()

    def get_node(self, p):
        for n in self.g.nodes():
            if points_are_close(p, n):
                return n
        return []

    def add_edge(self, n1, n2, entity):
        self.g.add_edge(n1, n2, object=entity)

    def get_edge(self, eg):
        return [[e[0], e[1], e[2]['object']] for e in self.g.edges(data=True)
                if e[2]['object'] is eg]

    def _remove_edge(self, n1, n2):
        self.g.remove_edge(n1, n2)

    def remove_edge(self, edge):
        e = self.get_edge(edge)
        assert(len(e) == 1)
        self.g.remove_edge(e[0][0], e[0][1])

    def remove_edges(self, edges):
        for e in edges:
            self.remove_edge(e)

    def add_line(self, n1, n2):
        line = Line(Element(start=n1, end=n2))
        add_or_join(self,
                    line.node1(ndec),
                    line.node2(ndec),
                    line,
                    self.rtol,
                    self.atol)

    def elements(self, type):
        """return lists of objects"""
        return [e[2]['object'] for e in self.g.edges(data=True)
                if isinstance(e[2]['object'], type)]

    def arcs(self):
        """return lists of arcs"""
        return self.elements(Arc)

    def lines(self):
        """return lists of lines"""
        return self.elements(Line)

    def circles(self):
        """return list of circle nodes"""
        return [n[1]['object'] for n in self.g.nodes(data=True)
                if n[1] and isinstance(n[1]['object'], Circle)]

    def virtual_nodes(self):
        nodes = []
        for e in self.elements(Shape):
            nodes += e.get_nodes()
        return nodes

    def angle_nodes(self, center, angle, rtol, atol):
        nodes = []
        for n in self.g.nodes():
            if points_are_close(center, n, rtol, atol):
                # Da gibt es keinen brauchbaren Winkel
                nodes.append(n)
            elif np.isclose(angle, alpha_line(center, n), rtol, atol):
                nodes.append(n)
        return nodes

    def radius_nodes(self, center, radius, rtol, atol):
        nodes = []
        for n in self.g.nodes():
            d = distance(center, n)
            if np.isclose(d, radius, rtol, atol):
                # Da gibt es keinen brauchbaren Winkel
                nodes.append(n)
        return nodes

    def get_corner_list(self, center, angle, rtol=1e-04, atol=1e-04):
        # Die Funktion liefert eine sortierte Liste aller Nodes auf einer
        # Linie als Corner-Objekte.
        corners = [Corner(center, c)
                   for c in self.angle_nodes(center, angle, rtol, atol)]
        if len(corners) == 1:
            logger.debug('the center is a corner')
            corners.append(Corner(center, tuple(center)))
        if len(corners) > 1:
            corners.sort()
        return corners

    def start_min_corner(self, i):
        return self.start_corners[0][i]

    def start_max_corner(self, i):
        return self.start_corners[-1][i]

    def repair_hull_line(self, center, angle):
        # We need to set our own tolerance range
        # to find the right points
        rtol = 1e-4
        atol = 1e-4

        corners = self.get_corner_list(center, angle, rtol, atol)
        if len(corners) < 2:
            # no hull without more than 1 corners
            logger.info('repair_hull_line: only {} corners'.
                        format(len(corners)))
            return

        [c.set_keep_node() for c in corners if c.is_equal(center, rtol, atol)]
        for p1, p2, data in [e for e in self.g.edges(data=True)]:
            clist_p1 = [c for c in corners if c.is_equal(p1, 0.0, 1e-7)]
            clist_p2 = [c for c in corners if c.is_equal(p2, 0.0, 1e-7)]
            if len(clist_p1) > 1:
                logger.warning("WARNING: {} corners and p1 close together"
                               .format(len(clist_p1)))
            if len(clist_p2) > 1:
                logger.warning("WARNING: {} corners and p2 close together"
                               .format(len(clist_p2)))

            if clist_p1 and clist_p2:
                # Both points are in the hull
                el = data['object']
                if isinstance(el, Line):
                    self._remove_edge(p1, p2)
                else:
                    [corner.set_keep_node() for corner in clist_p1]
                    [corner.set_keep_node() for corner in clist_p2]
                    if isinstance(el, Arc):
                        alpha_start = el.startangle
                        alpha_end = el.endangle
                        alpha_mid = middle_angle(alpha_start, alpha_end)
                        self._remove_edge(p1, p2)
                        a1 = Arc(Element(center=el.center,
                                         radius=el.radius,
                                         start_angle=alpha_start*180/np.pi,
                                         end_angle=alpha_mid*180/np.pi))
                        a2 = Arc(Element(center=el.center,
                                         radius=el.radius,
                                         start_angle=alpha_mid*180/np.pi,
                                         end_angle=alpha_end*180/np.pi))
                        self.add_edge(p1, a1.node2(ndec), a1)
                        self.add_edge(a2.node1(ndec), p2, a2)
            else:
                clist = []
                if clist_p1:
                    clist = [c for c in clist_p1]
                elif clist_p2:
                    clist = [c for c in clist_p2]
                [corner.set_keep_node() for corner in clist]

        [self.g.remove_node(c.point()) for c in corners if not c.keep_node()]

        # Rebuild Corner-list after correction
        corners = self.get_corner_list(center, angle, rtol, atol)

        if len(corners) > 1:
            corners.sort()
            p1 = corners[0].point()
            for c in corners[1:]:
                p2 = c.point()
                self.add_edge(p1, p2, Line(Element(start=p1, end=p2)))
                p1 = p2

    def set_minmax_radius(self, center):
        self.min_radius = 99999.0
        self.max_radius = 0.0
        for e in self.elements(Shape):
            mm_dist = e.minmax_from_center(center)
            self.min_radius = min(self.min_radius, mm_dist[0])
            self.max_radius = max(self.max_radius, mm_dist[1])

    def complete_hull_line(self, center, angle):
        if self.max_radius == 0.0:
            self.set_minmax_radius(center)

        corners = self.get_corner_list(center, angle)
        assert(corners)
        c_min = Corner(center, point(center, self.min_radius, angle, ndec))
        c_max = Corner(center, point(center, self.max_radius, angle, ndec))

        c_first = corners[0]
        if not c_min.is_same_corner(c_first):
            c_min.is_new_point = True
            p1 = c_min.point()
            p2 = c_first.point()
            self.add_edge(p1, p2, Line(Element(start=p1, end=p2)))

        c_last = corners[len(corners)-1]
        if not c_max.is_same_corner(c_last):
            c_max.is_new_point = True
            p2 = c_max.point()
            p1 = c_last.point()
            self.add_edge(p1, p2, Line(Element(start=p1, end=p2)))

        return (c_min, c_max)

    def complete_hull_arc(self, center, startangle, startcorner,
                          endangle, endcorner, radius):
        nodes = self.radius_nodes(center, radius, 1e-04, 1e-04)

        if startcorner.is_new_point:
            start_p = startcorner.point()
            nodes_sorted = [(distance(start_p, n), n) for n in nodes
                            if not points_are_close(start_p, n)]
            nodes_sorted.sort()
            p = nodes_sorted[0][1]
            angle_p = alpha_line(center, p)
            self.add_edge(start_p, p, Arc(
                Element(center=center, radius=radius,
                        start_angle=startangle*180/np.pi,
                        end_angle=angle_p*180/np.pi)))

        if endcorner.is_new_point:
            end_p = endcorner.point()
            nodes_sorted = [(distance(end_p, n), n) for n in nodes
                            if not points_are_close(end_p, n)]
            inx = len(nodes_sorted)-1
            p = nodes_sorted[inx][1]
            angle_p = alpha_line(center, p)
            self.add_edge(p, end_p, Arc(
                Element(center=center, radius=radius,
                        start_angle=angle_p*180/np.pi,
                        end_angle=endangle*180/np.pi)))

    def get_corner_nodes(self, center, angle):
        rtol = 1e-4
        atol = 1e-4

        corners = self.get_corner_list(center, angle, rtol, atol)
        if len(corners) < 2:
            return ()  # not enough corners
        return (corners[0].point(), corners[len(corners)-1].point())

    def point_lefthand_side(self, p1, p2):
        alpha = alpha_line(p2, p1)

        nbrs = [n for n in self.g.neighbors(p2)
                if not (points_are_close(n, p1) or points_are_close(n, p2))]
        if len(nbrs) == 0:
            return None  # unexpected end

        angles = []
        for p in nbrs:
            e_dict = self.g.get_edge_data(p2, p)
            assert(e_dict)
            e = e_dict['object']
            px = e.center_of_connection(ndec)
            alphax = alpha_line(p2, px)
            if np.isclose(alpha, alphax, 1e-05, 0.00001):
                # print("   >>> alpha={}, alphax={}".format(alpha, alphax))
                angles.append((0.0, p))
            else:
                # we don't move more than 180 degrees
                angles.append((alpha_angle(alpha, alphax), p))

        if len(angles) == 0:
            return None

        angles.sort()
        return angles[len(angles)-1][1]

    def get_new_area(self, start_p1, start_p2, solo):
        logger.debug('==> start of get_new_area({}, {})'
                     .format(start_p1, start_p2))
        e_dict = self.g.get_edge_data(start_p1, start_p2)
        if not e_dict:
            raise ValueError("Fatal: no edge-data found")
        area = []
        e = e_dict['object']
        x = e.get_point_number(start_p1)

        if e_dict[x]:
            logger.debug("<== area already tracked ({}) ***".format(x))
            return None
        e_dict[x] = True  # footprint
        area.append(e)

        if (isinstance(e, Circle) and not isinstance(e, Arc)):
            logger.debug("<== area is a circle ***")
            e_dict[1] = True  # footprint
            e_dict[2] = True  # footprint
            return area

        first_p = start_p1
        this_p = start_p2

        next_p = self.point_lefthand_side(first_p, this_p)
        if not next_p:
            logger.debug("<== dead end ({}, {})"
                         .format(first_p, this_p))
            return None

        a = normalise_angle(alpha_points(first_p, this_p, next_p))
        alpha = a

        afternext_p = ()
        c = 0
        while not (points_are_close(next_p, start_p1) and
                   points_are_close(afternext_p, start_p2)):
            logger.debug('  Next {}'.format(next_p))
            c += 1
            if c > self.num_edges * 2:
                logger.info("FATAL: *** over {} elements in area ? ***"
                            .format(self.num_edges))
                plot_area(area)
                sys.exit(1)
            e_dict = self.g.get_edge_data(this_p, next_p)
            e = e_dict['object']
            x = e.get_point_number(this_p)
            if e_dict[x]:
                logger.debug('<== path already tracked ***')
                return None
            e_dict[x] = True  # footprint
            first_p = this_p
            this_p = next_p
            next_p = self.point_lefthand_side(first_p, this_p)
            if not next_p:
                logger.debug("<== dead end ({},{})"
                             .format(first_p, this_p))
                return None

            a = normalise_angle(alpha_points(first_p, this_p, next_p))
            alpha += a
            area.append(e)
            afternext_p = self.point_lefthand_side(this_p, next_p)

        logger.debug("  END OF get_new_area")

        e_dict = self.g.get_edge_data(this_p, next_p)
        e = e_dict['object']
        x = e.get_point_number(this_p)
        e_dict[x] = True  # footprint
        area.append(e)
        a = normalise_angle(alpha_points(this_p, next_p, start_p2))
        alpha += a

        if alpha < 0.0:
            logger.debug("<== turn left expected, but it turned right ({})"
                         .format(alpha))
            return None

        logger.debug("<== area found !! ")
        return area

    def create_list_of_areas(self):
        """ return list of areas for each node and their neighbors
        """
        if len(self.area_list) > 0:
            # list already available
            return

        def append(area_list, a):
            for area in area_list:
                if area.is_identical(a):
                    return
            area_list.append(a)

        logger.debug("create new area list")

        if nxversion == 1:
            nx.set_edge_attributes(self.g, 0, True)
            nx.set_edge_attributes(self.g, 1, False)
            nx.set_edge_attributes(self.g, 2, False)
        else:
            nx.set_edge_attributes(self.g, True, 0)
            nx.set_edge_attributes(self.g, False, 1)
            nx.set_edge_attributes(self.g, False, 2)

        for p in self.g.nodes():
            if self.debug:
                print('.', end='', flush=True)
            neighbors = [n for n in self.g[p]]
            for next_p in neighbors:
                area = self.get_new_area(p, next_p, len(neighbors) < 3)
                if area:
                    a = Area(area, self.center, 0.0)
                    append(self.area_list, a)

        if self.debug:
            print(" done. {} areas found".format(len(self.area_list)))

    def list_of_areas(self):
        self.create_list_of_areas()
        return self.area_list

    def remove_all_areas(self):
        self.create_list_of_areas()
        for area in self.area_list:
            area.remove_edges(self.g, ndec)

    def copy_line(self, center, radius, start_angle, end_angle,
                  start_line, end_line, inner_circle, outer_circle, e,
                  rtol=1e-04,
                  atol=1e-04):
        """ Die Funktion kopiert die Teile einer Linie, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Line))
        if is_same_angle(start_angle, end_angle):
            points = inner_circle.intersect_line(e,
                                                 rtol,
                                                 atol,
                                                 False) + \
                     outer_circle.intersect_line(e,
                                                 rtol,
                                                 atol,
                                                 False) + \
                     [e.p2]
        else:
            points_start = e.intersect_line(start_line,
                                            rtol,
                                            atol,
                                            False)
            points_end = e.intersect_line(end_line,
                                          rtol,
                                          atol,
                                          False)
            points_inner = inner_circle.intersect_line(e,
                                                       rtol,
                                                       atol,
                                                       False)
            points_outer = outer_circle.intersect_line(e,
                                                       rtol,
                                                       atol,
                                                       False)
            points = points_start + points_end + \
                points_inner + points_outer + [e.p2]

        new_elements = []
        sorted_points = []

        for p in points:
            dist = distance(e.p1, p)
            sorted_points.append((dist, p))
        sorted_points.sort()

        p1 = e.p1
        for x, p2 in sorted_points:
            pm = middle_point_of_line(p1, p2)
            if is_point_inside_region(pm, center,
                                      inner_circle.radius,
                                      outer_circle.radius,
                                      start_angle, end_angle):
                new_elements.append(Line(Element(start=p1, end=p2)))
            p1 = p2

        return new_elements

    def copy_arc(self, center, radius, start_angle, end_angle,
                 start_line, end_line, inner_circle, outer_circle, e,
                 rtol=1e-04,
                 atol=1e-04):
        """ Die Funktion kopiert die Teile eines Kreissegments, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Arc))
        if is_same_angle(start_angle, end_angle):
            points = inner_circle.intersect_arc(e,
                                                rtol,
                                                atol,
                                                False) + \
                     outer_circle.intersect_arc(e,
                                                rtol,
                                                atol,
                                                False) + \
                     [e.p2]
        else:
            points_start = e.intersect_line(start_line,
                                            rtol,
                                            atol,
                                            False)
            points_end = e.intersect_line(end_line,
                                          rtol,
                                          atol,
                                          False)
            points_inner = inner_circle.intersect_arc(e,
                                                      rtol,
                                                      atol,
                                                      False)
            points_outer = outer_circle.intersect_arc(e,
                                                      rtol,
                                                      atol,
                                                      False)
            points = points_start + points_end + \
                points_inner + points_outer + [e.p2]

        new_elements = []
        sorted_points = []

        alpha_start = alpha_line(e.center, e.p1)

        for p in points:
            alpha_next = alpha_line(e.center, p)
            if less_equal(alpha_next, alpha_start):
                alpha_next += 2*np.pi
            sorted_points.append((alpha_next, p))
            alpha_start = alpha_next
        sorted_points.sort()

        p1 = e.p1
        alpha_start = alpha_line(e.center, e.p1)
        for x, p2 in sorted_points:
            alpha_end = alpha_line(e.center, p2)
            pm = middle_point_of_arc(e.center, e.radius, p1, p2, rtol=rtol)
            if is_point_inside_region(pm, center,
                                      inner_circle.radius, outer_circle.radius,
                                      start_angle, end_angle):

                if not (len(points) > 1 and
                        points_are_close(p1, p2, 1e-3, 1e-3)):
                    new_elements.append(
                        Arc(Element(center=e.center,
                                    radius=e.radius,
                                    start_angle=alpha_start*180/np.pi,
                                    end_angle=alpha_end*180/np.pi)))

            alpha_start = alpha_end
            p1 = p2
        return new_elements

    def copy_circle(self, center, radius, start_angle, end_angle,
                    start_line, end_line, inner_circle, outer_circle, e,
                    rtol=1e-04,
                    atol=1e-04):
        """ Die Funktion kopiert die Teile eines Kreises, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Circle))
        if is_same_angle(start_angle, end_angle):
            points = inner_circle.intersect_circle(e,
                                                   rtol,
                                                   atol,
                                                   False) + \
                     outer_circle.intersect_circle(e,
                                                   rtol,
                                                   atol,
                                                   False)
        else:
            points = e.intersect_line(start_line,
                                      rtol,
                                      atol) + \
                     e.intersect_line(end_line,
                                      rtol,
                                      atol) + \
                     inner_circle.intersect_circle(e,
                                                   rtol,
                                                   atol,
                                                   False) + \
                     outer_circle.intersect_circle(e,
                                                   rtol,
                                                   atol,
                                                   False)

        new_elements = []
        if len(points) < 2:
            if is_point_inside_region(e.p1, center,
                                      inner_circle.radius,
                                      outer_circle.radius,
                                      start_angle, end_angle):
                new_elements.append(
                    Circle(Element(center=e.center, radius=e.radius)))
            return new_elements

        sorted_points = []
        for p in points:
            alpha_p = alpha_line(e.center, p)
            sorted_points.append((alpha_p, p))
        sorted_points.sort()

        x, px = sorted_points[0]
        del sorted_points[0]
        p1 = px
        alpha_start = alpha_line(e.center, p1)
        for x, p2 in sorted_points:
            alpha_end = alpha_line(e.center, p2)
            pm = middle_point_of_arc(e.center, e.radius, p1, p2, rtol=rtol)
            if is_point_inside_region(pm, center,
                                      inner_circle.radius,
                                      outer_circle.radius,
                                      start_angle, end_angle):
                alpha_middle = middle_angle(alpha_start, alpha_end)
                arc1 = Arc(Element(center=e.center,
                                   radius=e.radius,
                                   start_angle=alpha_start*180/np.pi,
                                   end_angle=alpha_middle*180/np.pi))
                arc2 = Arc(Element(center=e.center,
                                   radius=e.radius,
                                   start_angle=alpha_middle*180/np.pi,
                                   end_angle=alpha_end*180/np.pi))
                new_elements.append(arc1)
                new_elements.append(arc2)

            alpha_start = alpha_end
            p1 = p2

        alpha_end = alpha_line(e.center, px)
        pm = middle_point_of_arc(e.center, e.radius, p1, px, rtol=rtol)
        if is_point_inside_region(pm, center,
                                  inner_circle.radius, outer_circle.radius,
                                  start_angle, end_angle):
            alpha_middle = middle_angle(alpha_start, alpha_end)
            arc1 = Arc(Element(center=e.center,
                               radius=e.radius,
                               start_angle=alpha_start*180/np.pi,
                               end_angle=alpha_middle*180/np.pi))
            arc2 = Arc(Element(center=e.center,
                               radius=e.radius,
                               start_angle=alpha_middle*180/np.pi,
                               end_angle=alpha_end*180/np.pi))
            new_elements.append(arc1)
            new_elements.append(arc2)
        return new_elements

    def copy_shape(self,
                   center,
                   radius,
                   startangle,
                   endangle,
                   inner_radius,
                   outer_radius,
                   split=False,
                   rtol=0.0,
                   atol=0.0):
        """ Die Funktion kopiert die Teile von Shape-Objekten, welche sich in
            der durch die Parameter definierten Teilkreisfläche befinden.
        """
        logger.debug('copy_shape({}, {})'.format(startangle, endangle))
        if not rtol:
            rtol = self.rtol
        if not atol:
            atol = self.atol

        if is_same_angle(startangle, endangle):
            start_line = Line(
                Element(start=center,
                        end=point(center, radius+1, startangle)))
            end_line = Line(
                Element(start=center,
                        end=point(center, radius+1, startangle)))
        else:
            start_line = Line(
                Element(start=center,
                        end=point(center, radius+1, startangle)))
            end_line = Line(
                Element(start=center,
                        end=point(center, radius+1, endangle)))

        if np.isclose(normalise_angle(startangle),
                      normalise_angle(endangle), 0.0):
            inner_circle = Circle(Element(center=center, radius=inner_radius))
            outer_circle = Circle(Element(center=center, radius=outer_radius))
        else:
            inner_circle = Arc(Element(center=center, radius=inner_radius,
                                       start_angle=startangle*180/np.pi,
                                       end_angle=endangle*180/np.pi))
            outer_circle = Arc(Element(center=center, radius=outer_radius,
                                       start_angle=startangle*180/np.pi,
                                       end_angle=endangle*180/np.pi))

        new_elements = []
        for e in self.elements(Shape):
            if isinstance(e, Line):
                new_elements += self.copy_line(
                    center, radius, startangle, endangle,
                    start_line, end_line,
                    inner_circle, outer_circle, e,
                    rtol=rtol,
                    atol=atol)

            elif isinstance(e, Arc):
                new_elements += self.copy_arc(
                    center, radius, startangle, endangle,
                    start_line, end_line,
                    inner_circle, outer_circle, e,
                    rtol=rtol,
                    atol=atol)

            elif isinstance(e, Circle):
                new_elements += self.copy_circle(
                    center, radius, startangle, endangle,
                    start_line, end_line,
                    inner_circle, outer_circle, e,
                    rtol=rtol,
                    atol=atol)

        if split:
            logger.debug('new Geometry with split')
            return Geometry(new_elements, 0.05, 0.1, split=split)
        else:
            return Geometry(new_elements, self.rtol, self.atol)

    def is_new_angle(self, alpha_list, alpha):
        for a in alpha_list:
            if np.isclose(a, alpha):
                return False
        return True

    def find_symmetry(self, center, radius,
                      startangle, endangle, sym_tolerance):
        logger.info("find symmetry")
        arealist = self.list_of_areas()

        logger.info("  {} areas available".format(len(arealist)))
        if len(arealist) == 0:
            return False

        arealist.sort()

        def add(areas, a):
            if not a.is_circle():
                for area in areas:
                    if area.is_equal(a, sym_tolerance):
                        area.increment(a)
                        return
            areas.append(a)

        arealist_match = []
        for a in arealist:
            a.sym_tolerance = sym_tolerance
            add(arealist_match, a)

        for a in arealist_match:
            a.set_delta()
        arealist_match.sort()

        area = arealist_match[0]
        if area.delta == 0.0:
            logger.info("No symmetry-axis found (delta == 0.0)")
            return False

        sym = part_of_circle(0.0, area.delta, 1)
        if sym == 0.0:
            logger.info("No symmetry-axis found (sym = 0.0)")
            return False
        area.delta = 2*np.pi/sym

        for alpha in area.symmetry_lines(startangle, endangle):
            p = point(center, radius+5, alpha)
            line = Line(Element(start=center, end=p))
            self.add_cut_line(line)

        self.sym_area = area
        return True

    def has_symmetry_area(self):
        return self.sym_area is not None

    def symmetry_startangle(self):
        return self.sym_area.sym_startangle

    def symmetry_endangle(self):
        return self.sym_area.sym_endangle

    def get_symmetry_copies(self):
        if self.sym_counterpart == 0:
            return int(self.sym_part)

        x = gcd(self.sym_part, self.sym_counterpart)
        return int(self.sym_part / x - 1)

    def is_mirrored(self):
        return len(self.mirror_corners) > 0

    def get_alfa(self):
        if self.is_mirrored():
            return 2*self.alfa
        else:
            return self.alfa

    def __str__(self):
        return "name...........: {}\n".format(self._name) + \
               "kind...........: {}\n".format(self.kind) + \
               "sym_part.......: {}\n".format(self.sym_part) + \
               "sym_counterpart: {}\n".format(self.sym_counterpart) + \
               "alpha..........: {}\n".format(self.alfa) + \
               "radius.........: {} -- {}\n".format(self.min_radius,
                                                    self.max_radius)

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self.__dict__ == other.__dict__)

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __hash__(self):
        """Override the default hash behavior
        (that returns the id or the object)"""
        return hash(tuple(sorted(self.__dict__.items())))

    def set_name(self, name):
        self._name = name

    def get_name(self):
        return self._name

    def minmax(self):
        """ find min and max coordinates of all shapes
            ([<min-x>, <max-x>, <min-y>, <max-y>])
        """
        mm = [99999, -99999, 99999, -99999]

        for e in self.elements(Shape):
            n = e.minmax()
            mm[0] = min(mm[0], n[0])
            mm[1] = max(mm[1], n[1])
            mm[2] = min(mm[2], n[2])
            mm[3] = max(mm[3], n[3])

        return mm

    def add_cut_line(self, e):
        self.cut_lines.append(e)

    def clear_cut_lines(self):
        self.cut_lines = []

    def render_cut_lines(self, renderer):
        for e in self.cut_lines:
            e.render(renderer, 'darkred')

    def render_airgaps(self, renderer):
        for g in self.airgaps:
            if isinstance(g, Arc):
                renderer.arc(g.startangle, g.endangle,
                             g.center, g.radius,
                             color='red')
            elif isinstance(g, Circle):
                renderer.circle(g.center, g.radius,
                                color='red')
            elif isinstance(g, Line):
                renderer.line(g.p1, g.p2,
                              color='red')
            elif isinstance(g, Point):
                renderer.point(g.p1, 'ro')

    def render_neighbors(self, renderer):
        for n in self.g.nodes():
            nbr_list = [nbr for nbr in self.g.neighbors(n)]
            if len(nbr_list) == 1:
                renderer.point(n, 'ro', color='orange')
            elif len(nbr_list) == 2:
                renderer.point(n, 'ro', color='green')
            elif len(nbr_list) == 3:
                renderer.point(n, 'ro', color='red')
            elif len(nbr_list) == 4:
                renderer.point(n, 'ro', color='blue')
            elif len(nbr_list) > 4:
                renderer.point(n, 'ro', color='black')

    def render_area_fill(self, renderer):
        legend = {}
        for area in self.list_of_areas():
            if area.is_iron():
                area.render_fill(renderer, 0.3)
                if area.name() and area.name() not in legend:
                    legend[area.name()] = area.render_legend(renderer, 0.3)

        for area in self.list_of_areas():
            if area.is_air():
                area.render_fill(renderer, 1.0)

        # magnet has no air inside
        for area in self.list_of_areas():
            if area.is_magnet():
                area.render_fill(renderer)
                if area.name() and area.name() not in legend:
                    legend[area.name()] = area.render_legend(renderer, 1.0)

        # winding has no air inside
        for area in self.list_of_areas():
            if area.is_winding():
                area.render_fill(renderer)
                if area.name() and area.name() not in legend:
                    legend[area.name()] = area.render_legend(renderer, 1.0)

        if legend:
            return [h for (k, h) in legend.items()]
        return []

    def get_points_in_iron(self):
        points = []
        for area in self.list_of_areas():
            p = area.get_point_inside(self)
            if p:
                points.append(p)
        return points

    def check_hull(self, center, radius, x, y, rtol, atol):
        node_count = 0
        miss_count = 0
        for h in convex_hull(self.virtual_nodes()):
            dist = distance(center, h)
            node_count += 1

            if not np.isclose(dist, radius, rtol, atol):
                if x is not None:
                    if np.isclose(x, h[0], rtol, atol):
                        continue
                if y is not None:
                    if np.isclose(y, h[1], rtol, atol):
                        continue
                miss_count += 1
        return miss_count == 0

    def get_machine(self):
        mm = self.minmax()
        height = mm[3]-mm[2]
        width = mm[1]-mm[0]

        c = []
        r = 0.0
        atol = 3.0

        if np.isclose(height, width, self.rtol, self.atol):
            r = width/2
            c = [mm[1]-r, mm[3]-r]
            logger.info("check for full machine")
            if self.check_hull(c, r, None, None, self.rtol, atol):
                logger.info("  it is full")
                return Machine(self, c, r, 0.0, 0.0)

            logger.info("check for quarter machine")
            r = width
            c = [mm[0], mm[2]]
            if self.check_hull(c, r, mm[0], mm[2], self.rtol, atol):
                logger.info("  it is a quarter")
                return Machine(self, c, r, 0.0, np.pi/2)

        elif np.isclose(width, height*2, self.rtol, self.atol):
            r = width/2
            c = [mm[1]-height, mm[2]]
            logger.info("check for half machine")
            if self.check_hull(c, r, None, mm[2], self.rtol, atol):
                logger.info("  it is a half")
                return Machine(self, c, r, 0.0, np.pi)

            c = [mm[1]-height, mm[3]]
            if self.check_hull(c, r, None, mm[3], self.rtol, atol):
                logger.info("  it is a half")
                return Machine(self, c, r, np.pi, 0.0)

        elif np.isclose(width*2, height, self.rtol, self.atol):
            r = width
            c = [mm[0], mm[1]-width]
            logger.info("check for half machine")
            c = [mm[1], mm[3]-width]
            if self.check_hull(c, r, mm[1], None, self.rtol, atol):
                logger.info("  it is a half")
                return Machine(self, c, r, np.pi/2.0, -np.pi/2.0)

            c = [mm[0], mm[3]-width]
            if self.check_hull(c, r, mm[0], None, self.rtol, atol):
                logger.info("  it is a half")
                return Machine(self, c, r, -np.pi/2.0, np.pi/2.0)

        machine = self.get_machine_part(mm)
        if machine:
            return machine

        logger.info("The shape of the Machine is unexpected")
        return Machine(self, [0.0, 0.0], 0.0, 0.0, 0.0)

    def is_new_center(self, center_list, center, rtol, atol):
        for c in center_list:
            if points_are_close(c[0], center, rtol, atol):
                c[1][0] += 1
                return False
        return True

    def get_machine_part(self, mm):
        center_list = []
        for e in self.elements(Arc):
            center = [round(e.center[0], 3), round(e.center[1], 3)]
            if self.is_new_center(center_list, center, self.rtol, self.atol):
                center_list.append((center, [1]))

        center = []
        count = 0
        unique = False
        for c in center_list:
            if c[1][0] == count:
                unique = False
            elif c[1][0] > count:
                count = c[1][0]
                unique = True
                center = c[0]

        if not unique:
            # Wir finden keine Arc-Objekte, welche uns einen Hinweis auf den
            # Center geben können. Wir versuchen in der Verzweiflung mit
            # x(min) und y(min)
            center = [round(mm[0], 4), round(mm[2], 4)]

        logger.debug("Center is {}".format(center))

        min_radius = 99999
        max_radius = 0
        startangle = 999.0
        endangle = -999.0

        for h in convex_hull(self.virtual_nodes()):
            angle = alpha_line(center, [round(h[0], 4), round(h[1], 4)])
            if angle < 0.0:
                logger.debug("strange point {}".format(h))
            startangle = min(startangle, angle)
            endangle = max(endangle, angle)

            dist = distance(center, h)
            min_radius = min(min_radius, dist)
            max_radius = max(max_radius, dist)

        center_left = round(center[0] - mm[0], 4)
        center_right = round(mm[1] - center[0], 4)
        center_down = round(center[1] - mm[2], 4)
        center_up = round(mm[3] - center[1], 4)

        min_r = min(center_left, center_right, center_up, center_down)
        max_r = max(center_left, center_right, center_up, center_down)

        if less_equal(center_left, 0.0) or less_equal(center_right, 0.0) or \
           less_equal(center_up, 0.0) or less_equal(center_down, 0.0):
            if not np.isclose(center_left, center_right):
                x = part_of_circle(startangle, endangle)

                if x > 2:
                    return Machine(self,
                                   [round(center[0], 8), round(center[1], 8)],
                                   max_radius, startangle, endangle)

        if min_radius >= max_radius*0.9 and min_r >= max_r*0.9:
            # Mit 10 % Abweichungen gehen wir noch von einem ganzen Motor aus.
            return Machine(self,
                           [round(center[0], 8), round(center[1], 8)],
                           max(max_radius, max_r), 0.0, 0.0)

        if np.isclose(center_down, 0.0):
            min_r = min(center_left, center_right, center_up)
            if min_r >= max_r*0.9:
                # Vermutlich ein halber Motor
                return Machine(self,
                               [round(center[0], 8), round(center[1], 8)],
                               max(max_radius, max_r), 0.0, np.pi)

            if np.isclose(center_left, 0.0):
                min_r = min(center_right, center_up)
                if min_r >= max_r*0.9:
                    # Vermutlich ein viertel Motor
                    return Machine(self,
                                   [round(center[0], 8), round(center[1], 8)],
                                   max(max_radius, max_r), 0.0, np.pi/2)

        # TODO: handle half and quarter machines

        if np.isclose(center_left, center_right):
            # Der x-Wert des Center scheint zu stimmen.
            # Eine Symmetrie an der y-Achse ist möglich.
            # Dem y-Wert ist nicht zu trauen!
            if center_down > 0.0:
                # This center is too high.
                nodes = [n for n in convex_hull(
                    self.g.nodes()) if n[0] > center[0]]
                assert(nodes)
                p = nodes[0]
                for n in nodes[1:]:
                    if n[1] < p[1]:
                        p = n

                m_min = 99999.0
                for n in nodes:
                    m = line_m(p, n)
                    if m and m > 0.0:
                        # print("m = {}, p={}, n={}".format(m, p, n))
                        m_min = min(m_min, m)

                y = line_n([p[0]-center[0], p[1]], m_min)
                center[1] = y
                angle = alpha_line(center, p)

        return Machine(self, [round(center[0], 8), round(center[1], 8)],
                       max_radius, angle, np.pi - angle)

    def is_new_radius(self, radius_list, radius):
        for r, d in radius_list:
            if np.isclose(r, radius):
                return False
        return True

    def is_airgap(self, center, radius, startangle, endangle, circle, atol):
        """ Die Funktion untersucht, ob sich der Parameter circle in einem
            Luftspalt befindet.
        """
        ok = True
        borders = 0
        for e in self.elements(Shape):
            for p in e.intersect_circle(circle, 0.0, atol, True):
                if not self.is_border_line(center,
                                           startangle, endangle,
                                           e, atol):
                    logger.info("BAD: Point {}".format(p))
                    self.airgaps.append(Point(p))
                    self.airgaps.append(e)
                    ok = False
                else:
                    borders += 1
        return (ok, borders)

    def is_border_line(self, center, startangle, endangle, e, atol):
        if isinstance(e, Line):
            if np.isclose(startangle, endangle):
                return False  # full

            angle_p1 = alpha_line(center, e.p1)
            if np.isclose(startangle, angle_p1, 1e-3, atol):
                angle_p2 = alpha_line(center, e.p2)
                return np.isclose(startangle, angle_p2, 1e-3, atol)
            elif np.isclose(endangle, angle_p1, 1e-3, atol):
                angle_p2 = alpha_line(center, e.p2)
                return np.isclose(endangle, angle_p2, 1e-3, atol)
        return False

    def get_gaplist(self, center):
        gaplist = []
        for e in self.elements(Shape):
            gaplist += [e.minmax_from_center(center)]
        gaplist.sort()

        airgaps = []
        dist_max = 0.0
        for g in gaplist:
            if not less_equal(g[0], dist_max):
                airgaps.append((dist_max, g[0]))
            dist_max = max(dist_max, g[1])
        return airgaps

    def detect_airgaps(self, center, startangle, endangle, atol):
        """ Die Funktion sucht Luftspalt-Kandidaten und liefert eine Liste
            von Möglichkeiten mit jeweils einem minimalen und einem maximalen
            Radius als Begrenzung des Luftspalts.
        """
        gaplist = []
        for e in self.elements(Shape):
            if not self.is_border_line(center, startangle, endangle, e, atol):
                gaplist += [e.minmax_from_center(center)]

        gaplist.sort()

        self.min_radius = gaplist[0][0]
        self.max_radius = gaplist[-1][1]

        airgaps = []
        min_radius = self.min_radius + 1.0
        cur_radius = gaplist[0][1]
        max_radius = self.max_radius - 1.0

        for g in gaplist:
            if greater(g[0], cur_radius) and \
               greater(cur_radius, min_radius) and \
               less(g[0], max_radius):
                airgaps.append((cur_radius, g[0]))

            cur_radius = max(cur_radius, g[1])

        return airgaps

    def get_circles(self, center, radius):
        return [c for c in self.elements(Circle)
                if points_are_close(center, c.center) and
                np.isclose(radius, c.radius)]

    def delete_circle(self, center, radius):
        for c in self.elements(Circle):
            if points_are_close(center, c.center):
                if np.isclose(radius, c.radius):
                    logger.info("Center circle removed")
                    self.remove_edge(c)
                    return

    def alpha_of_circles(self, circles, center):
        angle = 0.0
        for c in circles:
            if isinstance(c, Arc):
                alpha_c_p1 = alpha_line(center, c.p1)
                alpha_c_p2 = alpha_line(center, c.p2)
                angle += alpha_angle(alpha_c_p1, alpha_c_p2)
            else:
                angle = 2*np.pi
        return angle

    def delete_airgap_circle(self, center,
                             lower_radius, radius, upper_radius,
                             angle_tot):
        lower_circles = self.get_circles(center, lower_radius)
        angle_sum = self.alpha_of_circles(lower_circles, center)
        if angle_sum / angle_tot < 0.75:
            return False

        upper_circles = self.get_circles(center, upper_radius)
        angle_sum = self.alpha_of_circles(upper_circles, center)
        if angle_sum / angle_tot < 0.75:
            return False

        self.remove_edges(self.get_circles(center, radius))
        return True

    def create_auxiliary_lines(self, rightangle, leftangle):
        for area in self.list_of_areas():
            logger.debug("create_auxiliary_lines for {}"
                         .format(area.identifier()))

            areas_inside = [a for a in self.area_list
                            if area.is_inside(a)]
            areas_border = {a.identifier(): a for a in areas_inside
                            if area.has_connection(self, a, ndec)}
            areas_notouch = {a.identifier(): a for a in areas_inside
                             if not area.has_connection(self, a, ndec)}

            for id, a in areas_notouch.items():
                for id2, a2 in areas_border.items():
                    if a.has_connection(self, a2, ndec):
                        areas_border[id2] = a2
                        break

            for id, a in areas_border.items():
                if id in areas_notouch:
                    del areas_notouch[id]

            notouch_list = []
            for id, a in areas_notouch.items():
                if not notouch_list:
                    notouch_list.append({id: a})
                else:
                    touched = False
                    for l in notouch_list:
                        for id2, a2 in l.items():
                            if a.has_connection(self, a2, ndec):
                                l[id] = a
                                touched = True
                                break
                        if touched:
                            break
                    if not touched:
                        notouch_list.append({id: a})

            for lst in notouch_list:
                if lst:
                    gap_list = []
                    for id, a in lst.items():
                        logger.debug(">> area {} not in touch"
                                     .format(a.identifier()))
                        gap_list += area.get_lowest_gap_list(a,
                                                             self.center,
                                                             self.max_radius,
                                                             rightangle,
                                                             leftangle)
                    gap_list.sort()
                    assert(len(gap_list) > 0)
                    points = gap_list[0][1]

                    logger.debug(">> auxiliary line from {} to {}"
                                 .format(points[0], points[1]))
                    line = Line(Element(start=points[0], end=points[1]),
                                color='orange',
                                attr='auxline')
                    add_or_join(self, points[0], points[1], line,
                                self.rtol, self.atol)

    def set_rotor(self):
        self.sym_counterpart = 1
        self.sym_part = 0

    def set_stator(self):
        self.sym_counterpart = 1
        self.sym_part = 2

    def is_rotor(self):
        if self.sym_counterpart:
            return self.sym_part < self.sym_counterpart
        return False

    def is_stator(self):
        if self.sym_counterpart:
            return self.sym_part > self.sym_counterpart
        return False

    def num_variable(self):
        if self.is_stator():
            return 'm.num_sl_gen'
        if self.is_rotor():
            return 'm.npols_gen'
        return 'm.{}_ncopies'.format(self.kind)

    def _delete_a_tiny_element(self, n0, n1, dict01, n2):
        dict12 = self.g.get_edge_data(n1, n2)
        if dict12.get('deleted', False):
            return False

        n12_el = dict12['object']
        if not isinstance(n12_el, Line):
            return False

        logger.debug("tiny line from {} to {} deleted"
                     .format(n0, n2))
        self._remove_edge(n0, n1)
        self._remove_edge(n1, n2)
        dict12['deleted'] = True
        dict01['deleted'] = True
        line = Line(Element(start=n0, end=n2))
        self.add_edge(n0, n2, line)
        return True

    def delete_tiny_elements(self, mindist):
        if mindist == 0.0:
            return

        edges = [edge for edge in self.g.edges(data=True)
                 if distance(edge[0], edge[1]) < mindist]
        if len(edges) > 0:
            logger.info("mindist={}: {} tiny elements found"
                        .format(mindist, len(edges)))
        deleted = 0
        for edge in edges:
            if edge[2].get('deleted', False):
                continue

            nbrs_n1 = [nbr for nbr in self.g.neighbors(edge[0])
                       if nbr is not edge[1]]
            nbrs_n2 = [nbr for nbr in self.g.neighbors(edge[1])
                       if nbr is not edge[0]]

            if len(nbrs_n1) == 1:
                if self._delete_a_tiny_element(edge[1], edge[0],
                                               edge[2], nbrs_n1[0]):
                    deleted += 1
                    continue

            if len(nbrs_n2) == 1:
                if self._delete_a_tiny_element(edge[0], edge[1],
                                               edge[2], nbrs_n2[0]):
                    deleted += 1
                    continue

        logger.info("{} tiny elements deleted".format(deleted))
        return

    def search_subregions(self):
        if self.is_stator():
            return self.search_stator_subregions()

        if self.is_rotor():
            return self.search_rotor_subregions()

        logger.warning("no stator or rotor assigned")
        return self.search_unknown_subregions()

    def search_stator_subregions(self, place=''):
        is_inner = self.is_inner
        if place == 'in':
            is_inner = True
        elif place == 'out':
            is_inner = False

        for area in self.list_of_areas():
            area.mark_stator_subregions(is_inner,
                                        self.is_mirrored(),
                                        self.alfa,
                                        self.center,
                                        self.min_radius,
                                        self.max_radius)

        # windings close to endangle?
        wdg_areas = [a for a in self.list_of_areas()
                     if a.type == 2 and a.close_to_endangle]

        # air or iron near windings?
        air_areas = [a for a in self.list_of_areas() if a.type == 9]
        for a in air_areas:
            if a.around_windings(self.list_of_areas()):
                a.type = 0  # air
            elif a.close_to_endangle:
                if wdg_areas:
                    a.type = 0  # air
                else:
                    a.type = 6  # iron shaft (Zahn)
            else:
                a.type = 6  # iron shaft (Zahn)

        # yoke or shaft ?
        iron_areas = [a for a in self.list_of_areas() if a.type == 5]
        for a in iron_areas:
            if a.around_windings(self.list_of_areas()):
                a.type = 6  # iron shaft (Zahn)

    def search_rotor_subregions(self, place=''):
        is_inner = self.is_inner
        if place == 'in':
            is_inner = True
        elif place == 'out':
            is_inner = False

        types = {}
        for area in self.list_of_areas():
            t = area.mark_rotor_subregions(is_inner,
                                           self.is_mirrored(),
                                           self.alfa,
                                           self.center,
                                           self.min_radius,
                                           self.max_radius)
            if t in types:
                types[t] += 1
            else:
                types[t] = 1

        if 4 in types:  # magnet rectangle
            if types[4] > 1:
                logger.debug("{} magnets in rotor".format(types[4]))

            for area in self.list_of_areas():
                if area.type == 3:
                    area.type = 1  # iron

    def search_unknown_subregions(self):
        for area in self.list_of_areas():
            area.mark_unknown_subregions(self.is_mirrored(),
                                         self.alfa,
                                         self.center,
                                         self.min_radius,
                                         self.max_radius)
        return

    def num_areas_of_type(self, type):
        return len([area for area in self.list_of_areas()
                    if area.type == type])

    def num_of_windings(self):
        return self.num_areas_of_type(2)

    def area_close_to_endangle(self, type):
        return len([area for area in self.list_of_areas()
                    if area.type == type and area.close_to_endangle])

    def print_nodes(self):
        print("=== List of Nodes ({}) ===".format(self.number_of_nodes()))
        for n in self.g.nodes():
            print(n)

    def write_nodes(self, filename):
        nodes = []
        for n in self.g.nodes():
            nodes.append(n)
        nodes.sort()

        content = []
        for n in nodes:
            content.append(u"{}".format(n))
            for nbr in self.g.neighbors(n):
                content.append(u" --> {}".format(nbr))

        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write('\n'.join(content))
