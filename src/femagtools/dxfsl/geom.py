# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.geom
  ~~~~~~~~~~~~~~~~~~~~~

  a geometry is composed of connected shapes (lines, arcs, circles) which
  build areas

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function

import numpy as np
import networkx as nx
import logging
import inspect
import sys
from pathlib import Path
from femagtools.dxfsl.corner import Corner
from femagtools.dxfsl.area import Area
import femagtools.dxfsl.area as AREA
from .shape import Element, Shape, Circle, Arc, Line, Point
from .shape import is_Circle, is_Arc, is_Line
from .machine import Machine
from femagtools.dxfsl.concat import Concatenation
from femagtools.dxfsl.areabuilder import AreaBuilder
from femagtools.dxfsl.functions import Timer
from femagtools.dxfsl.journal import Journal, getJournal
from .functions import less_equal, less, greater, greater_equal
from .functions import distance, alpha_line, alpha_points, alpha_angle
from .functions import point, points_are_close, is_point_inside_region
from .functions import line_m, line_n, lines_intersect_point
from .functions import middle_point_of_line, middle_point_of_arc
from .functions import middle_angle, positive_angle
from .functions import normalise_angle, is_same_angle
from .functions import part_of_circle, gcd
from .functions import point_on_arc, points_on_line, nodes_are_equal
from .functions import area_size
import io
import time

logger = logging.getLogger('femagtools.geom')


nxversion = int(nx.__version__.split('.')[0])
geom_mindist = 0.0


def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno


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

def create_geometry(new_elements, split=False):
    return Geometry(new_elements,
                    center=(0.0, 0.0),
                    split=split)


def intersect_and_split(inp_elements, rtol, atol):
    logger.info("Load input elements ... %s", len(inp_elements))
    out_elements = []
    for e in inp_elements:
        out_size = len(out_elements)
        intersect_and_split_element(e, out_elements, 0, out_size, rtol, atol)
    logger.debug(" done", e)
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
    split_el = el.overlapping_shape(out_elements[x], rtol, atol)
    if split_el:
        logger.debug("=== overlapping elements ===")
        out_elements[x] = None
        return split_el

    points = el.intersect_shape(out_elements[x], rtol, atol, True)
    if len(points) > 0:
        split_elements = out_elements[x].split(points, rtol, atol)
        if len(split_elements) > 0:
            out_elements += split_elements
            out_elements[x] = None
        split_el = el.split(points, rtol, atol)
        return split_el
    return []


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

nodes_filecount = 0

TYPE_UNDEFINED = 0
TYPE_ROTOR = 1
TYPE_STATOR = 2

class Geometry(object):
    """collection of connected shapes"""

    def __init__(self,
                 elements=[],
                 center=[],
                 rtol=1e-04,
                 atol=1e-03,
                 is_inner=False,
                 is_outer=False,
                 split=False,
                 concatenate=False,
                 connect=False,
                 delete=False,
                 adjust=False,
                 main=False,
                 type=TYPE_UNDEFINED,
                 debug=False):
        self._name = ''
        self.kind = ''
        self.mirror_corners = []
        self.start_corners = []
        self.end_corners = []
        self.sym_part = 0
        self.sym_counterpart = 0
        self.sym_type = type
        self.sym_slices = 0
        self.sym_slice_angle = 0.0
        self.alfa = 0.0
        self.center = []
        self.with_center_node = False
        self.min_radius = 0.0
        self.max_radius = 0.0
        self.is_inner = is_inner
        self.is_outer = is_outer
        self.cut_lines = []
        self.sym_area = None
        self.airgaps = []
        self.area_list = []
        self.areagroup_list = []
        self.g = nx.Graph()
        self.rtol = rtol
        self.atol = atol
        self.debug = debug
        self.num_edges = 0
        self.wdg_is_mirrored = False
        self.has_windings = False
        self.has_magnets = False
        self.journal = getJournal()
        self.area_errors = 0
        self.critical_points = []
        i = 0

        logger.debug("Geometry(split=%s, concat=%s, connect=%s, delete=%s, adjust=%s, main=%s,",
                     split, concatenate, connect, delete, adjust, main)
        logger.debug("         rtol=%s, atol=%s)",
                     rtol, atol)

        timer = Timer(start_it=True)
        self.c_concat = 0
        self.c_connect = 0

        def get_elements(elements, split):
            if split:
                return intersect_and_split(elements, self.rtol, self.atol)

            src_elements = [e for e in elements]
            if main:
                self.journal.put_elements(len(src_elements))

            if not concatenate:
                return src_elements

            concat = Concatenation(rtol=self.rtol, atol=self.atol)
            c, new_elements = concat.concatenate_matching_elements(src_elements,
                                                                   main=main)
            self.c_concat = c
            return new_elements

        nbr_nodes = []
        if concatenate:
            elements = [e for e in elements]
            geom = Geometry(elements, center=(0.0, 0.0))
            logger.debug("REAL GEOMETRY START ")
            elements = geom.copy_all_elements(alpha=0.0)
            nbr_nodes = geom.get_nodes(num_of_nbrs=[2, 4])

        for e in get_elements(elements, split):
            if e:
                e.id = i
                n = self.find_nodes(e.start(), e.end())
                try:
                    self.add_or_join_edge(n[0], n[1], e,
                                          rtol=self.rtol,
                                          atol=self.atol)
                except Exception as ex:
                    logger.warn("EXCEPTION %s", ex)
                    if e:  # must be a circle
                        self.g.add_node(e.center, object=e)
            i += 1

        self.num_edges = self.number_of_edges()

        if center:
            self.set_center(center)

        if connect:
            self.c_connect = self.connect_all_appendices(rtol=1e-04,
                                                         atol=1e-03,
                                                         main=True)
            if delete:
                self.delete_all_appendices()

            if concatenate and self.c_concat > 0:
                self.c_connect += self.connect_all_nodes(rtol=1e-04,
                                                         atol=1e-03,
                                                         additional_nodes=nbr_nodes,
                                                         main=True)

        if adjust:
            self.adjust_all_points()

        if delete and center:
            self.set_minmax_radius()

        timer.stop("-- Geometry initialised in %0.4f seconds --")
        logger.debug("End Geometry(concatenated=%s, connected=%s)",
                     self.c_concat, self.c_connect)

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
            e[2]['object'].move(offset, ndec)
        for c in self.circles():
            c.move((-offset[0], -offset[1]))
        mapping = {n: (round(n[0] + offset[0], ndec),
                       round(n[1] + offset[1], ndec))
                   for n in self.g.nodes()}
        nx.relabel_nodes(self.g, mapping, copy=False)

    def rotate_nodes(self, alpha, nodes):
        T = np.array(((np.cos(alpha), -np.sin(alpha)),
                      (np.sin(alpha), np.cos(alpha))))
        rotnodes = np.dot(T, np.asarray(nodes).T).T.tolist()
        return rotnodes

    def rotate(self, alpha):
        """rotates all objects by angle alpha"""
        logger.debug("rotate geometry(%s)", alpha)
        T = np.array(((np.cos(alpha), -np.sin(alpha)),
                      (np.sin(alpha), np.cos(alpha))))
        for e in self.g.edges(data=True):
            e[2]['object'].transform(T, alpha, ndec)
        for c in self.circles():
            c.transform(T, alpha, ndec)
        rotnodes = np.dot(T, np.asarray(self.g.nodes()).T).T.tolist()
        mapping = {n: (round(r[0], ndec),
                       round(r[1], ndec))
                   for n, r in zip(self.g.nodes(), rotnodes)}

        nx.relabel_nodes(self.g, mapping, copy=False)

    def check_geom(self, what):
        logger.debug("check geometry of %s", what)
        parts = int(round(np.pi * 2 / self.alfa, 1))
        logger.debug(" --parts......: %s", parts)
        logger.debug(" --ist alpha..: %s", self.alfa)
        real_alfa = np.pi * 2 / parts
        logger.debug(" --soll alpha.: %s", real_alfa)
        if not np.isclose(self.alfa, real_alfa):
            logger.debug(" --BAD angle ==> get corrected machine")
            return self.correct_geom(real_alfa)
        return None

    def correct_geom(self, correct_alpha):
        elements = []
        for e in self.g.edges(data=True):
            o = e[2]['object'].correct(self.alfa, correct_alpha, ndec)
            elements.append(o)
        geom = Geometry(elements,
                        center=self.center,
                        is_inner=self.is_inner,
                        is_outer=self.is_outer,
                        rtol=self.rtol,
                        atol=self.atol,
                        type=self.sym_type)
        geom.alfa = correct_alpha
        geom.kind = self.kind
        geom.sym_part = self.sym_part
        geom.sym_type = self.sym_type
        return geom

    def log_geom(self):
        logger.info("Kind............: %s", self.kind)
        logger.info("Center..........: %s", self.center)
        logger.info("Alpha...........: %s", self.alfa)
        logger.info("Is Inner........: %s", self.is_inner)
        logger.info("Is Outer........: %s", self.is_outer)
        logger.info("Min Radius......: %s", self.min_radius)
        logger.info("Max Radius......: %s", self.max_radius)
        logger.info("Mirror Corners..: %s", self.mirror_corners)
        logger.info("Start Corners...: %s", self.start_corners)
        logger.info("End Corners.....: %s", self.end_corners)
        logger.info("Edges...........: %s", self.num_edges)
        if self.sym_slices > 0:
            logger.info("Symmetry Slices.: %s", self.sym_slices)
            logger.info("Slice Angle.....: %s", self.sym_slice_angle)

    def scale(self, factor):
        """scales all objects"""
        for e in self.g.edges(data=True):
            e[2]['object'].scale(factor)
        for c in self.circles():
            c.scale(factor)
        mapping = {n: (round(factor * n[0], ndec),
                       round(factor * n[1], ndec))
                   for n in self.g.nodes()}
        nx.relabel_nodes(self.g, mapping, copy=False)
        self.diameters = tuple([factor*d for d in self.diameters])

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

    def find_other_node(self, node, **kwargs):
        """return closest node to node in arg within pickdist"""
        nodes = list(kwargs.get('g', self.g))
        if nodes:
            anodes = np.asarray(nodes)
            # la.norm on numpy below 1.8 does not accept axis
            c = anodes - node
            dist = np.sqrt(np.einsum('ij, ij->i', c, c))
            idx = dist.argmin()
            candidates = [(dist[i], nodes[i]) for i in range(len(dist))
                          if dist[i] < 0.01]
            if not candidates:
                return None
            candidates.sort()
            if len(candidates) > 1:
                i = 0
                d, n = candidates[0]
                if d == 0.0:
                    d, n = candidates[1]
                return n
        return None

    def find_the_node(self, p, **kwargs):
        """return closest nodes to points in arg within pickdist"""
        nodes = list(kwargs.get('g', self.g))
        if nodes:
            anodes = np.asarray(nodes)
            # la.norm on numpy below 1.8 does not accept axis
            c = anodes - p
            dist = np.sqrt(np.einsum('ij, ij->i', c, c))
            # dist = la.norm(np.asarray(nodes) - p, axis=1)
            idx = dist.argmin()
            if dist[idx] < self.atol:
                return nodes[idx]
        return None

    def polyline_paths(self):
        """return lists of line paths"""
        paths = []
        g = self.g.copy()
        g.remove_edges_from(self.arcs())
        while self.number_of_edges() > 0:
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

    def add_element(self, e, rtol, atol):
        n = self.find_nodes(e.start(), e.end())
        try:
            self.add_or_join_edge(n[0], n[1], e, rtol=rtol, atol=atol)
        except Exception as ex:
            logger.warn("EXCEPTION in add_element: %s", ex)

    def add_or_join_edge(self, n1, n2, entity, rtol=1e-03, atol=1e-03):
        """ adds a new entity to graph or joins entity with existing
        geom: Geometry
        n1, n2: nodes
        entity
        """
        if n1 == n2:
            logger.debug(
                "add_or_join_edge: Tiny element with same node on both sides ignored: %s", n1)
            logger.debug(
                "add_or_join_edge: -- element: %s", entity)
            return

        e = self.get_edge_element(n1, n2)
        if not e:  # no duplicates
            self.add_edge(n1, n2, entity)
            return

        logger.debug("add_or_join_edge: Duplicate connection: %s <--> %s", n1, n2)
        if is_Line(e):
            if is_Line(entity):
                logger.debug("add_or_join_edge: Duplicate Lines ignored")
                return  # its ok

        if is_Arc(e):
            if is_Arc(entity):
                if points_are_close(e.center, entity.center, rtol=rtol, atol=atol):
                    if points_are_close(e.p1, entity.p1):
                        logger.debug("add_or_join_edge: Duplicate Arcs ignored")
                        return  # its ok

        if is_Circle(entity):
            if is_Circle(e):
                logger.debug("add_or_join_edge: Duplicate Circle ignored")
                return  # its ok

        if is_Circle(entity) or is_Circle(e):
            e1, e2 = entity.cut_into_halves()
            logger.debug("add_or_join_edge: Element near circle is cut into halves")
            self.add_element(e1, rtol, atol)
            self.add_element(e2, rtol, atol)
            return  # halves installed

        m1 = e.center_of_connection()
        m2 = entity.center_of_connection()
        logger.debug("add_or_join_edge: midpoints: %s -- %s", m1, m2)
        if points_are_close(m1, m2, rtol, 1e-2):
            logger.debug("add_or_join_edge: Elements are close together")
            return  # ok

        e1, e2 = entity.cut_into_halves()
        logger.debug("add_or_join_edge: cut into halves")
        self.add_element(e1, rtol, atol)
        self.add_element(e2, rtol, atol)
        return  # halves installed

    def add_edge(self, n1, n2, entity):
        if points_are_close(n1, n2):
            logger.debug("WARNING in add_edge(): Points of %s are close together",
                         entity.classname())
            logger.debug("        n1 = %s, n2 = %s", n1, n2)
            logger.debug("        p1 = %s, p2 = %s", entity.p1, entity.p2)

        if self.has_edge(n1, n2):
            logger.warning("FATAL ERROR: Duplicates in add_edge(%s, %s)", n1, n2)

        entity.set_nodes(n1, n2)
        logger.debug("add_edge %s - %s  (%s)", n1, n2, entity.classname())
        self.g.add_edge(n1, n2, object=entity)

    def has_edge(self, n1, n2):
        if self.g.get_edge_data(n1, n2):
            return True
        return False

    def get_edge(self, obj):
        return [[e[0], e[1], e[2]['object']] for e in self.g.edges(data=True)
                if e[2]['object'] is obj]

    def get_edge_element(self, n1, n2):
        e_dict = self.g.get_edge_data(n1, n2)
        if e_dict:
            return e_dict.get('object', None)
        return None

    def get_edge_nodes(self, edge):
        e = self.get_edge(edge)
        assert(len(e) == 1)
        return (e[0][0], e[0][1])

    def _remove_edge(self, n1, n2):
        logger.debug("remove_edge %s - %s", n1, n2)
        self.g.remove_edge(n1, n2)
        self._remove_node(n1)
        self._remove_node(n2)

    def remove_edge(self, edge):
        e = self.get_edge(edge)
        if len(e) != 1:
            logger.info("remove edge failed: %s", edge)
            raise ValueError("remove edge failed")
        assert(len(e) == 1)
        self._remove_edge(e[0][0], e[0][1])

    def remove_edges(self, objs):
        for o in objs:
            self.remove_edge(o)

    def _remove_node(self, n):
        for nbr in self.g.neighbors(n):
            return
        try:
            self.g.remove_node(n)
        except nx.NetworkXError:
            logger.warning("WARNING: remove node %s failed", n)

    def add_line(self, n1, n2, color=None, linestyle=None):
        line = Line(Element(start=n1, end=n2),
                    color=color,
                    linestyle=linestyle)
        self.add_element(line,
                         rtol=self.rtol,
                         atol=self.atol)

    def add_arc(self, n1, n2, center, radius, color=None, linestyle=None):
        angle_n1 = alpha_line(center, n1)
        angle_n2 = alpha_line(center, n2)
        arc = Arc(Element(center=center,
                          radius=radius,
                          start_angle=angle_n1*180/np.pi,
                          end_angle=angle_n2*180/np.pi),
                  color=color,
                  linestyle=linestyle)
        self.add_element(arc,
                         rtol=self.rtol,
                         atol=self.atol)

    def elements(self, type=Shape):
        """return lists of objects"""
        return [e[2]['object'] for e in self.g.edges(data=True)
                if isinstance(e[2]['object'], type)]

    def elements_and_nodes(self, type=Shape):
        """return lists of objects"""
        return [(e[0], e[1], e[2]['object']) for e in self.g.edges(data=True)
                if isinstance(e[2]['object'], type)]

    def arcs(self):
        """return lists of arcs"""
        return self.elements(Arc)

    def lines(self):
        """return lists of lines"""
        return self.elements(Line)

    def split_lines_longer_than(self, length):
        """split lines longer than length"""
        new_lines = []
        rem_lines = []
        for p1, p2, data in [e for e in self.g.edges(data=True)
                             if isinstance(e[2]['object'], Line)]:
            ln = data['object']
            if ln.length() > length:
                p = ln.center_of_connection()
                new_lines += ln.split([p])
                rem_lines.append((p1, p2))

        for p1, p2 in rem_lines:
            self._remove_edge(p1, p2)
        for new_ln in new_lines:
            self.add_element(new_ln,
                             rtol=self.rtol,
                             atol=self.atol)

    def split_all_lines_longer_than(self, length):
        """split lines longer than length"""
        new_lines = []
        rem_lines = []
        elist = self.elements_and_nodes(Line)

        for p1, p2, line in elist:
            line_len = line.length()
            if line_len < length:
                continue
            d = abs(distance(self.center, p1) - distance(self.center, p2))
            parts = 3
            if d > line_len / 2:
                parts = 5
            elif d > line_len / 3:
                parts = 4
            p_start = p1

            for p_next in points_on_line(p1, p2, parts=parts):
                new_lines.append(Line(Element(start=p_start, end=p_next)))
                p_start = p_next

            new_lines.append(Line(Element(start=p_start, end=p2)))
            rem_lines.append((p1, p2))

        for p1, p2 in rem_lines:
            self._remove_edge(p1, p2)
        for new_line in new_lines:
            self.add_element(new_line,
                             rtol=self.rtol,
                             atol=self.atol)

    def circles(self):
        """return list of circle nodes"""
        return [n[1]['object'] for n in self.g.nodes(data=True)
                if n[1] and isinstance(n[1]['object'], Circle)]

    def nodes(self):
        for n in self.g.nodes():
            yield n

    def get_nodes(self, num_of_nbrs=[]):
        if num_of_nbrs:
            nodes = []
            for n in self.g.nodes():
                nbr_list = [nbr for nbr in self.g.neighbors(n)]
                if len(nbr_list) in num_of_nbrs:
                    nodes.append(n)
        else:
            nodes = [n for n in self.g.nodes()]
        return nodes

    def virtual_nodes(self):
        nodes = []
        for e in self.elements(Shape):
            nodes += e.get_nodes()
        return nodes

    def get_neighbors(self, n):
        return [nbr for nbr in self.g.neighbors(n)]

    def num_of_neighbors(self, n):
        nbrs = [nbr for nbr in self.g.neighbors(n)]
        return len(nbrs)

    def angle_nodes(self, center, angle, rtol, atol):
        if np.isclose(abs(angle), np.pi, rtol, atol):
            angle_func = positive_angle
        else:
            angle_func = normalise_angle
        angle = angle_func(angle)

        nodes = []
        for n in self.g.nodes():
            if points_are_close(center, n, rtol, atol):
                # Da gibt es keinen brauchbaren Winkel
                nodes.append(n)
            else:
                angle_line = angle_func(alpha_line(center, n))
                if np.isclose(angle, angle_line, rtol, atol):
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
        center_added = len(corners) == 1
        if len(corners) == 1:
            logger.debug('get_corner_list: the center is a corner')
            corners.append(Corner(center, tuple(center)))
        if len(corners) > 1:
            corners.sort()
        return center_added, corners

    def start_min_corner(self, i):
        return self.start_corners[0][i]

    def start_max_corner(self, i):
        return self.start_corners[-1][i]

    def get_start_airgap_node(self):
        if self.is_inner:
            return self.start_corners[-1]
        else:
            return self.start_corners[0]

    def dist_start_max_corner(self):
        logger.debug("begin of dist_start_max_corner")
        logger.debug("start corners: %s", self.start_corners)
        d = distance(self.center, self.start_corners[-1])
        logger.debug("end of dist_start_max_corner: %s", d)
        return d

    def dist_end_max_corner(self, mirrored=True):
        logger.debug("begin of dist_end_max_corner")
        logger.debug("end corners: %s", self.end_corners)

        if self.is_mirrored() and mirrored:
            return self.dist_start_max_corner()
        d = distance(self.center, self.end_corners[-1])
        logger.debug("end of dist_end_max_corner: %s", d)
        return d

    def max_corners_match(self):
        d1 = self.dist_start_max_corner()
        d2 = self.dist_end_max_corner(mirrored=False)
        return np.isclose(d1, d2, rtol=1e-3, atol=1e-3)

    def dist_start_min_corner(self):
        logger.debug("begin of dist_start_min_corner")
        logger.debug("start corners: %s", self.start_corners)
        d = distance(self.center, self.start_corners[0])
        logger.debug("end of dist_start_min_corner: %s", d)
        return d

    def dist_end_min_corner(self, mirrored=True):
        logger.debug("begin of dist_end_min_corner")
        logger.debug("end corners: %s", self.end_corners)

        if self.is_mirrored() and mirrored:
            return self.dist_start_min_corner()
        d = distance(self.center, self.end_corners[0])
        logger.debug("end of dist_end_min_corner: %s", d)
        return d

    def min_corners_match(self):
        d1 = self.dist_start_min_corner()
        d2 = self.dist_end_min_corner(mirrored=False)
        return np.isclose(d1, d2, rtol=1e-3, atol=1e-3)

    def min_max_corners_match(self):
        return self.min_corners_match() and self.max_corners_match()

    def get_start_airgap_corner(self):
        if self.is_inner:
            p = (self.max_radius, 0.0)
            cp = self.start_corners[-1]
        else:
            p = (self.min_radius, 0.0)
            cp = self.start_corners[0]
        if points_are_close(p, cp, rtol=1e-2, atol=1e-1):
            if points_are_close(p, cp, rtol=1e-3, atol=1e-2):
                logger.debug("get_start_airgap_corner: critical")
                logger.debug(" -- soll: %s,  ist: %s", p, cp)
            return cp, True
        return p, False

    def get_end_airgap_corner(self):
        if self.is_inner:
            p = point(self.center, self.max_radius, self.alfa, ndec)
            cp = self.end_corners[-1]
        else:
            p = point(self.center, self.min_radius, self.alfa, ndec)
            cp = self.end_corners[0]
        logger.debug("End Airgap Corner: %s is %s", p, cp)
        if points_are_close(p, cp, rtol=1e-2, atol=1e-1):
            if points_are_close(p, cp, rtol=1e-3, atol=1e-2):
                logger.debug("get_end_airgap_corner: critical")
                logger.debug(" -- soll: %s,  ist: %s", p, cp)
            return cp, True
        return p, False

    def get_start_airgap_corner_point(self):
        p, b = self.get_start_airgap_corner()
        return Point(p)

    def get_end_airgap_corner_point(self):
        p, b = self.get_end_airgap_corner()
        return Point(p)

    def get_airgap_radius(self):
        if self.is_inner:
            return self.max_radius
        if self.is_outer:
            return self.min_radius
        return None

    def get_opposite_radius(self):
        if self.is_inner:
            return self.min_radius
        if self.is_outer:
            return self.max_radius
        return None

    def area_size(self):
        pts = [p for p in self.start_corners]
        end_pts = [p for p in reversed(self.end_corners)]
        return area_size(pts + end_pts)

    def repair_hull_line(self, center, angle, corners, with_center, rtol=None, atol=None):
        # We need to set our own tolerance range
        # to find the right points
        if not rtol:
            rtol = 1e-3
        if not atol:
            atol = 1e-3

        logger.debug("begin repair_hull_line(center=%s, angle=%s)", center, angle)
        [logger.debug(" --> Corner %s", c) for c in corners]

        if len(corners) < 2:
            # no hull without more than 1 corners
            logger.debug('end of repair_hull_line: only %s corners: EXIT',
                         len(corners))
            return

        [c.set_keep_node() for c in corners if c.is_equal(center, rtol, atol)]
        for p1, p2, data in [e for e in self.g.edges(data=True)]:
            clist_p1 = [c for c in corners if c.is_equal(p1, 0.0, 1e-7)]
            clist_p2 = [c for c in corners if c.is_equal(p2, 0.0, 1e-7)]
            if len(clist_p1) > 1:
                logger.warning("WARNING in repair_hull_line(): %s corners and p1 close together",
                               len(clist_p1))
                logger.warning("        p1 = %s, p2 = %s", p1, p2)
            if len(clist_p2) > 1:
                logger.warning("WARNING in repair_hull_line(): %s corners and p2 close together",
                               len(clist_p2))
                logger.warning("        p2 = %s, p1 = %s", p2, p1)

            if clist_p1 and clist_p2:
                # Both points are in the hull
                el = data['object']
                if isinstance(el, Line):
                    logger.debug("remove Line: %s <> %s", p1, p2)
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
                        self.add_element(a1,
                                         rtol=self.rtol,
                                         atol=self.atol)
                        self.add_element(a2,
                                         rtol=self.rtol,
                                         atol=self.atol)
            else:
                clist = []
                if clist_p1:
                    clist = [c for c in clist_p1]
                elif clist_p2:
                    clist = [c for c in clist_p2]
                [corner.set_keep_node() for corner in clist]

        # Rebuild Corner-list after correction
        center_added, corners = self.get_corner_list(center, angle, rtol, atol)
        for c in corners:
            logger.debug("Correct Corner: %s", c)

        if with_center or self.with_center_node:
            c_corner = Corner(center, tuple(center))
            if c_corner not in corners:
                corners.append(c_corner)

        if len(corners) > 1:
            corners.sort()
            p1 = corners[0].point()
            for c in corners[1:]:
                p2 = c.point()
                self.add_element(Line(Element(start=p1, end=p2)),
                                 rtol=self.rtol,
                                 atol=self.atol)
                p1 = p2

        self.set_minmax_radius()
        logger.debug('end of repair_hull_line')

    def set_center(self, center):
        self.center = center
        self.set_minmax_radius()

    def set_minmax_radius(self):
        if not self.center:
            raise ValueError("FATAL ERROR: no center in Geometry")
        self.min_radius = 99999.0
        self.max_radius = 0.0
        for e in self.elements(Shape):
            min_dist, max_dist = e.minmax_from_center(self.center)
            self.min_radius = min(self.min_radius, min_dist)
            self.max_radius = max(self.max_radius, max_dist)

    def complete_hull_line(self, angle):
        logger.debug("begin complete_hull_line")
        if not self.center:
            raise ValueError("FATAL ERROR: no center in Geometry")

        center_added, corners = self.get_corner_list(self.center, angle)
        assert(corners)
        c_min = Corner(self.center, point(self.center, self.min_radius, angle, ndec))
        c_max = Corner(self.center, point(self.center, self.max_radius, angle, ndec))

        c_first = corners[0]
        if not c_min.is_same_corner(c_first):
            c_min.is_new_point = True
            p1 = c_min.point()
            p2 = c_first.point()
            self.add_element(Line(Element(start=p1, end=p2)),
                             rtol=self.rtol,
                             atol=self.atol)

        c_last = corners[len(corners)-1]
        if not c_max.is_same_corner(c_last):
            c_max.is_new_point = True
            p2 = c_max.point()
            p1 = c_last.point()
            self.add_element(Line(Element(start=p1, end=p2)),
                             rtol=self.rtol,
                             atol=self.atol)
        logger.debug("end complete_hull_line")
        return (c_min, c_max)

    def complete_hull_arc(self, startangle, startcorner,
                          endangle, endcorner, radius):
        nodes = self.radius_nodes(self.center, radius, 1e-04, 1e-04)

        if startcorner.is_new_point:
            start_p = startcorner.point()
            nodes_sorted = [(distance(start_p, n), n) for n in nodes
                            if not points_are_close(start_p, n)]
            nodes_sorted.sort()
            p = nodes_sorted[0][1]
            angle_p = alpha_line(self.center, p)
            arc = Arc(Element(center=self.center, radius=radius,
                              start_angle=startangle*180/np.pi,
                              end_angle=angle_p*180/np.pi))
            self.add_element(arc,
                             rtol=self.rtol,
                             atol=self.atol)

        if endcorner.is_new_point:
            end_p = endcorner.point()
            nodes_sorted = [(distance(end_p, n), n) for n in nodes
                            if not points_are_close(end_p, n)]
            inx = len(nodes_sorted)-1
            p = nodes_sorted[inx][1]
            angle_p = alpha_line(self.center, p)
            arc = Arc(Element(center=self.center, radius=radius,
                              start_angle=angle_p*180/np.pi,
                              end_angle=endangle*180/np.pi))
            self.add_element(p, end_p, arc,
                             rtol=self.rtol,
                             atol=self.atol)

    def get_corner_nodes(self, center, angle):
        rtol = 1e-4
        atol = 1e-4

        center_added, corners = self.get_corner_list(center, angle, rtol, atol)
        if len(corners) < 2:
            return ()  # not enough corners
        return (corners[0].point(), corners[len(corners)-1].point())

    def set_start_corners(self, center, angle):
        self.start_corners = self.get_corner_nodes(center, angle)

    def set_end_corners(self, center, angle):
        self.end_corners = self.get_corner_nodes(center, angle)

    def get_angle(self, alpha1, alpha2):
        if np.isclose(alpha1, alpha2, 0.001, 0.001):
            return 0.0
        return alpha_angle(alpha1, alpha2)

    def set_edge_attributes(self):
        if nxversion == 1:
            nx.set_edge_attributes(self.g, 0, True)
            nx.set_edge_attributes(self.g, 1, False)
            nx.set_edge_attributes(self.g, 2, False)
        else:
            nx.set_edge_attributes(self.g, True, 0)
            nx.set_edge_attributes(self.g, False, 1)
            nx.set_edge_attributes(self.g, False, 2)

    def create_list_of_areas(self, main=False, delete=False):
        """ return list of areas for each node and their neighbors
        """
        if delete:  # clear list of areas
            self.area_list = []

        if len(self.area_list) > 0:
            logger.debug("area list already available")
            # list already available
            return

        areabuilder = AreaBuilder(geom=self)
        areabuilder.create_list_of_areas(main=main)
        self.area_list = areabuilder.area_list
        logger.debug("area list created")

    def list_of_areas(self):
        self.create_list_of_areas()
        return self.area_list

    def remove_all_areas(self):
        self.create_list_of_areas()
        for area in self.area_list:
            area.remove_edges(self.g, ndec)

    def intersect_the_line(self, line, rtol=1e-04, atol=1e-04):
        for e in self.elements(Shape):
            pts = e.intersect_line(line, rtol, atol, False)
            if pts:
                return pts
        return []

    def the_point_is_inside(self, p, rtol=1e-04, atol=1e-04):
        for e in self.elements(Shape):
            if e.is_point_inside(p, rtol=rtol, atol=atol, include_end=True):
                # logger.info("point %s is inside %s", p, e)
                return True
        return False

    def copy_line(self, center, radius, start_angle, end_angle,
                  start_line, end_line, inner_circle, outer_circle, e,
                  rtol=1e-04,
                  atol=1e-04,
                  points_inner=None,
                  points_outer=None):
        """ Die Funktion kopiert die Teile einer Linie, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Line))
        if is_same_angle(start_angle, end_angle):
            pts_inner = inner_circle.intersect_line(e,
                                                    rtol,
                                                    atol,
                                                    False)
            pts_outer = outer_circle.intersect_line(e,
                                                    rtol,
                                                    atol,
                                                    False)
            points = pts_inner + pts_outer + [e.p2]
        else:
            pts_start = e.intersect_line(start_line,
                                         rtol,
                                         atol,
                                         False)
            pts_end = e.intersect_line(end_line,
                                       rtol,
                                       atol,
                                       False)
            pts_inner = inner_circle.intersect_line(e,
                                                    rtol,
                                                    atol,
                                                    False)
            pts_outer = outer_circle.intersect_line(e,
                                                    rtol,
                                                    atol,
                                                    False)
            points = pts_start + pts_end + \
                pts_inner + pts_outer + [e.p2]

        if points_inner is not None and pts_inner:
            points_inner += pts_inner
        if points_outer is not None and pts_outer:
            points_outer += pts_outer

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
                 atol=1e-04,
                 points_inner=None,
                 points_outer=None):
        """ Die Funktion kopiert die Teile eines Kreissegments, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Arc))
        if is_same_angle(start_angle, end_angle):
            pts_inner = inner_circle.intersect_arc(e,
                                                   rtol,
                                                   atol,
                                                   False)
            pts_outer = outer_circle.intersect_arc(e,
                                                   rtol,
                                                   atol,
                                                   False)
            points = pts_inner + pts_outer + [e.p2]
        else:
            pts_start = e.intersect_line(start_line,
                                         rtol,
                                         atol,
                                         False)
            pts_end = e.intersect_line(end_line,
                                       rtol,
                                       atol,
                                       False)
            pts_inner = inner_circle.intersect_arc(e,
                                                   rtol,
                                                   atol,
                                                   False)
            pts_outer = outer_circle.intersect_arc(e,
                                                   rtol,
                                                   atol,
                                                   False)
            points = pts_start + pts_end + \
                pts_inner + pts_outer + [e.p2]

        if points_inner is not None and pts_inner:
            points_inner += pts_inner
        if points_outer is not None and pts_outer:
            points_outer += pts_outer

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
            pm = middle_point_of_arc(e.center, e.radius, p1, p2, rtol=rtol, atol=atol)
            if is_point_inside_region(pm, center,
                                      inner_circle.radius,
                                      outer_circle.radius,
                                      start_angle, end_angle,
                                      rtol=rtol, atol=atol):
                if not (len(points) > 1 and
                        points_are_close(p1, p2, rtol=rtol, atol=atol)):
                    if len(points) == 1 and e.rtheta is not None:
                        a = Arc(Element(center=e.center,
                                        radius=e.radius,
                                        start_angle=alpha_start*180/np.pi,
                                        end_angle=alpha_end*180/np.pi,
                                        width=e.width,
                                        height=e.height,
                                        rtheta=e.rtheta,
                                        start_param=e.start_param,
                                        end_param=e.end_param))
                    else:
                        a = Arc(Element(center=e.center,
                                        radius=e.radius,
                                        start_angle=alpha_start*180/np.pi,
                                        end_angle=alpha_end*180/np.pi))
                    if points_are_close(a.p1, a.p2, rtol=1e-02, atol=1e-02):
                        logger.debug("ATTENTION: creation of a tiny arc")
                        logger.debug("-- %s", a)
                        a.set_attribute("tiny")
                    if points_are_close(a.p1, a.p2, rtol=1e-06, atol=1e-06):
                        logger.debug("-- points are equal")
                    else:
                        new_elements.append(a)
            alpha_start = alpha_end
            p1 = p2
        return new_elements

    def copy_circle(self, center, radius, start_angle, end_angle,
                    start_line, end_line, inner_circle, outer_circle, e,
                    rtol=1e-04,
                    atol=1e-04,
                    points_inner=None,
                    points_outer=None):
        """ Die Funktion kopiert die Teile eines Kreises, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Circle))
        if is_same_angle(start_angle, end_angle):
            pts_inner = inner_circle.intersect_circle(e,
                                                      rtol,
                                                      atol,
                                                      False)
            pts_outer = outer_circle.intersect_circle(e,
                                                      rtol,
                                                      atol,
                                                      False)
            points = pts_inner + pts_outer
        else:
            pts_start = e.intersect_line(start_line,
                                         rtol,
                                         atol)
            pts_end = e.intersect_line(end_line,
                                       rtol,
                                       atol)
            pts_inner = inner_circle.intersect_circle(e,
                                                      rtol,
                                                      atol,
                                                      False)
            pts_outer = outer_circle.intersect_circle(e,
                                                      rtol,
                                                      atol,
                                                      False)
            points = pts_start + pts_end + pts_inner + pts_outer

        if points_inner is not None and pts_inner:
            points_inner += pts_inner
        if points_outer is not None and pts_outer:
            points_outer += pts_outer

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
        sorted_points.append((x, px))
        p1 = px
        alpha_start = alpha_line(e.center, p1)
        for x, p2 in sorted_points[1:]:
            alpha_end = alpha_line(e.center, p2)
            pm = middle_point_of_arc(e.center, e.radius, p1, p2, rtol=rtol, atol=atol)
            if is_point_inside_region(pm, center,
                                      inner_circle.radius,
                                      outer_circle.radius,
                                      start_angle, end_angle,
                                      rtol=rtol, atol=atol):
                a = Arc(Element(center=e.center,
                                radius=e.radius,
                                start_angle=alpha_start*180/np.pi,
                                end_angle=alpha_end*180/np.pi))
                new_elements.append(a)
            alpha_start = alpha_end
            p1 = p2
        return new_elements

    def copy_shape(self,
                   radius,
                   startangle,
                   endangle,
                   inner_radius,
                   outer_radius,
                   split=False,
                   rtol=0.0,
                   atol=0.0,
                   append_inner=False,
                   append_outer=False,
                   delete_appendices=False,
                   concatenate=True,
                   connect=True):
        """ Die Funktion kopiert die Teile von Shape-Objekten, welche sich in
            der durch die Parameter definierten Teilkreisfläche befinden.
        """
        logger.debug('begin copy_shape(%s, %s)', startangle, endangle)

        if not rtol:
            rtol = self.rtol
        if not atol:
            atol = self.atol
        rtol = 1e-5
        atol = 1e-4
        logger.debug(' -> rtol=%s,  atol=%s', rtol, atol)

        self.with_center_node = self.find_the_node(self.center) is not None

        if is_same_angle(startangle, endangle):
            start_line = Line(
                Element(start=self.center,
                        end=point(self.center, radius+1, startangle)))
            end_line = Line(
                Element(start=self.center,
                        end=point(self.center, radius+1, startangle)))
        else:
            start_line = Line(
                Element(start=self.center,
                        end=point(self.center, radius+1, startangle)))
            end_line = Line(
                Element(start=self.center,
                        end=point(self.center, radius+1, endangle)))

        if np.isclose(normalise_angle(startangle),
                      normalise_angle(endangle), 0.0):
            inner_circle = Circle(Element(center=self.center, radius=inner_radius))
            outer_circle = Circle(Element(center=self.center, radius=outer_radius))
        else:
            inner_circle = Arc(Element(center=self.center, radius=inner_radius,
                                       start_angle=startangle*180/np.pi,
                                       end_angle=endangle*180/np.pi))
            outer_circle = Arc(Element(center=self.center, radius=outer_radius,
                                       start_angle=startangle*180/np.pi,
                                       end_angle=endangle*180/np.pi))

        new_elements = []
        pts_inner = [] if append_inner else None
        pts_outer = [] if append_outer else None

        for e in self.elements(Shape):
            if isinstance(e, Line):
                new_elements += self.copy_line(
                    self.center, radius, startangle, endangle,
                    start_line, end_line,
                    inner_circle, outer_circle, e,
                    rtol=rtol,
                    atol=atol,
                    points_inner=pts_inner,
                    points_outer=pts_outer)

            elif isinstance(e, Arc):
                new_elements += self.copy_arc(
                    self.center, radius, startangle, endangle,
                    start_line, end_line,
                    inner_circle, outer_circle, e,
                    rtol=rtol,
                    atol=atol,
                    points_inner=pts_inner,
                    points_outer=pts_outer)
            elif isinstance(e, Circle):
                new_elements += self.copy_circle(
                    self.center, radius, startangle, endangle,
                    start_line, end_line,
                    inner_circle, outer_circle, e,
                    rtol=rtol,
                    atol=atol,
                    points_inner=pts_inner,
                    points_outer=pts_outer)

        if pts_inner and len(pts_inner) > 1:
            pts_inner.sort(reverse=True)
            p1 = pts_inner[0]
            for p2 in pts_inner[1:]:
                start_angle = alpha_line(self.center, p1)
                end_angle = alpha_line(self.center, p2)
                arc = Arc(Element(center=self.center,
                                  radius=inner_radius,
                                  start_angle=start_angle*180/np.pi,
                                  end_angle=end_angle*180/np.pi))
                new_elements.append(arc)
                p1 = p2

        if pts_outer and len(pts_outer) > 1:
            pts_outer.sort(reverse=True)
            p1 = pts_outer[0]
            for p2 in pts_outer[1:]:
                start_angle = alpha_line(self.center, p1)
                end_angle = alpha_line(self.center, p2)
                arc = Arc(Element(center=self.center,
                                  radius=outer_radius,
                                  start_angle=start_angle*180/np.pi,
                                  end_angle=end_angle*180/np.pi))
                new_elements.append(arc)
                p1 = p2

        center = self.center

        if split:
            logger.debug('new Geometry with split')

        geom = Geometry(new_elements,
                        center=center,
                        rtol=self.rtol,
                        atol=self.atol,
                        is_inner=self.is_inner,
                        is_outer=self.is_outer,
                        concatenate=concatenate,
                        connect=connect,
                        delete=delete_appendices,
                        split=split,
                        type=self.sym_type)
        geom.with_center_node = self.with_center_node

        logger.debug('end copy_shape')
        return geom

    def copy_all_elements(self, alpha):
        logger.debug("begin copy_all_elements(alpha=%s)", alpha)
        if alpha == 0.0:
            T = None
        else:
            T = np.array(((np.cos(alpha), -np.sin(alpha)),
                          (np.sin(alpha), np.cos(alpha))))

        all_el = []
        lines = 0
        arcs = 0
        circles = 0
        el = None
        for e in self.elements(Shape):
            if isinstance(e, Line):
                lines += 1
                el = Line(Element(start=e.p1,
                                  end=e.p2))
            elif isinstance(e, Arc):
                arcs += 1
                alpha_start = alpha_line(e.center, e.p1)
                alpha_end = alpha_line(e.center, e.p2)
                el = Arc(Element(center=e.center,
                                 radius=e.radius,
                                 start_angle=alpha_start*180/np.pi,
                                 end_angle=alpha_end*180/np.pi))
            elif isinstance(e, Circle):
                circles += 1
                el = Circle(Element(center=e.center,
                                    radius=e.radius))
            else:
                el = None
            if el is not None:
                if T is not None:
                    el.transform(T, alpha, ndec)
                all_el.append(el)

        logger.debug("end copy_all_elements: %s lines, %s arcs, %s circles",
                     lines, arcs, circles)
        return all_el

    def new_clone(self, new_elements,
                  split=False,
                  concatenate=False,
                  connect=False,
                  adjust=False):
        return Geometry(new_elements,
                        center=self.center,
                        rtol=self.rtol,
                        atol=self.atol,
                        is_inner=self.is_inner,
                        is_outer=self.is_outer,
                        split=split,
                        concatenate=concatenate,
                        connect=connect,
                        adjust=adjust,
                        type=self.sym_type)

    def is_new_angle(self, alpha_list, alpha):
        for a in alpha_list:
            if np.isclose(a, alpha):
                return False
        return True

    def find_symmetry(self, radius,
                      startangle, endangle, sym_tolerance):
        arealist = self.list_of_areas()
        logger.debug("begin of find_symmetry: - %s areas available", len(arealist))
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

        arealist_sym = [a for a in arealist_match if a.symmetry > 0]
        if not arealist_sym:
            logger.debug("end of find_symmetry: No symmetry-axis found (delta == 0.0)")
            return False

        ggt = arealist_sym[0].symmetry
        for a in arealist_sym[1:]:
            if ggt != a.symmetry:
                ggt = gcd(ggt, a.symmetry)
                if ggt == 1:
                    logger.warning("asymmetrical iteration of areas detected")
                    break

                if ggt != a.symmetry:
                    logger.warning("unhandled asymmetry")
                    break

        arealist_ok = [a for a in arealist_sym if a.symmetry == ggt]
        midlist = {}
        for a in arealist_ok:
            mid = round(a.start, 3)
            if midlist.get(mid, None) is None:
                midlist[mid] = [1, a]
            else:
                midlist[mid][0] = midlist[mid][0]+1

        arealist_srt = [[v[0], k, v[1]] for k, v in midlist.items()]
        arealist_srt.sort(reverse=True)

        if not arealist_srt:
            logger.debug("end of find_symmetry: no sorted arealist")
            return False

        area = arealist_srt[0][2]
        sym = area.symmetry
        area.delta = 2*np.pi/sym
        self.sym_slices = sym
        self.sym_slice_angle = area.delta

        for alpha in area.symmetry_lines(startangle, endangle):
            p = point(self.center, radius+5, alpha)
            line = Line(Element(start=self.center, end=p))
            self.add_cut_line(line)

        self.sym_area = area
        return True

    def has_symmetry_area(self):
        return self.sym_area is not None

    def symmetry_startangle(self):
        return self.sym_area.sym_startangle

    def symmetry_endangle(self):
        return self.sym_area.sym_endangle

    def rotate_symmetry_parameters(self):
        self.sym_startangle += self.sym_slice_angle
        self.sym_endangle += self.sym_slice_angle
        if self.sym_area:
            self.sym_area.sym_startangle = self.sym_startangle
            self.sym_area.sym_endangle = self.sym_endangle

    def get_symmetry_copies(self):
        if self.sym_counterpart == 0:
            return int(self.sym_part)

        x = gcd(self.sym_part, self.sym_counterpart)
        return int(self.sym_part / x - 1)

    def is_mirrored(self):
        return len(self.mirror_corners) > 0

    def winding_is_mirrored(self):
        return self.wdg_is_mirrored

    def get_alfa(self):
        if self.is_mirrored():
            return 2*self.alfa
        else:
            return self.alfa

    def __str__(self):
        real_alfa = 0
        if self.sym_part > 0:
            if self.is_mirrored():
                real_alfa = np.pi / self.sym_part
            else:
                real_alfa = 2*np.pi / self.sym_part
        return "name...........: {}\n".format(self._name) + \
               "kind...........: {}\n".format(self.kind) + \
               "sym_part.......: {}\n".format(self.sym_part) + \
               "sym_counterpart: {}\n".format(self.sym_counterpart) + \
               "alpha..........: {}\n".format(self.alfa) + \
               "alpha real.....: {}\n".format(real_alfa) + \
               "circle.........: {}\n".format(self.alfa * self.sym_part) + \
               "mirrored.......: {}\n".format(self.is_mirrored()) + \
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
                renderer.point(g.p1, 'o', color='red')

    def render_neighbors(self, renderer):
        for n in self.g.nodes():
            nbr_list = [nbr for nbr in self.g.neighbors(n)]
            if len(nbr_list) == 1:
                renderer.point(n, 'o', color='orange')
            elif len(nbr_list) == 2:
                renderer.point(n, 'o', color='green')
            elif len(nbr_list) == 3:
                renderer.point(n, 'o', color='red')
            elif len(nbr_list) == 4:
                renderer.point(n, 'o', color='blue')
            elif len(nbr_list) > 4:
                renderer.point(n, 'o', color='black')

    def render_area_fill(self, renderer):
        legend = {}
        for area in self.list_of_areas():
            if area.is_iron():
                area.render_fill(renderer)
                if area.name() and area.name() not in legend:
                    legend[area.name()] = area.render_legend(renderer)
            if area.is_shaft():
                area.render_fill(renderer)
                if area.name() and area.name() not in legend:
                    legend[area.name()] = area.render_legend(renderer)

        for area in self.list_of_areas():
            if area.is_air():
                area.render_fill(renderer)

        # magnet has no air inside
        for area in self.list_of_areas():
            if area.is_magnet():
                area.render_fill(renderer)
                if area.name() and area.name() not in legend:
                    legend[area.name()] = area.render_legend(renderer)

        # winding has no air inside
        for area in self.list_of_areas():
            if area.is_winding():
                area.render_fill(renderer)
                if area.name() and area.name() not in legend:
                    legend[area.name()] = area.render_legend(renderer)

        if legend:
            return [h for (k, h) in legend.items()]
        return []

    def reduce_winding_nodes(self, mindist=0.01):
        return self.reduce_element_nodes(mindist=mindist,
                                         area_types=(AREA.TYPE_WINDINGS,))

    def reduce_element_nodes(self, mindist=0.01, area_types=()):
        timer = Timer(start_it=True)
        nodes_deleted = 0
        for area in self.list_of_areas():
            if not area_types or area.type in area_types:
                nodes_deleted += area.reduce_element_nodes(self, mindist)

        t = timer.stop("-- {} nodes deleted in %0.4f seconds --".format(nodes_deleted))
        self.journal.put('time_deleting_nodes', t)
        if nodes_deleted:
            self.journal.put('nodes_deleted', nodes_deleted)
            self.area_list = []
        return nodes_deleted > 0

    def render_areagroups(self, renderer):
        if not self.areagroup_list:
            return
        for area in self.areagroup_list:
            area.render(renderer,
                        color="yellow",
                        fill=False)
        return

    def render_magnet_phi(self, renderer):
        magnets = [a for a in self.list_of_areas()]
        if not magnets:
            return
        arrow_len = [a.magnet_arrow_length() for a in magnets]
        length = max(arrow_len)

        for area in magnets:
            area.render_magnet_phi(renderer, length)

    def render_critical(self, renderer):
        for e in self.critical_points:
            e.render(renderer, 'darkred')

    def get_points_in_iron(self):
        points = []
        for area in self.list_of_areas():
            p = area.get_point_inside(self)
            if p:
                points.append(p)
        return points

    def check_hull(self, radius, x, y, rtol, atol):
        node_count = 0
        miss_count = 0
        for h in convex_hull(self.virtual_nodes()):
            dist = distance(self.center, h)
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
        atol = 3.0

        logger.debug("*** Begin of get_machine() ***")

        logger.debug(">> minmax: %s", mm)
        logger.debug(">> w=%s,  h=%s", width, height)

        if np.isclose(height, width, self.rtol, self.atol):
            radius = width/2
            self.set_center([mm[1]-radius, mm[3]-radius])
            logger.info("check for full machine")
            if self.check_hull(radius, None, None, self.rtol, atol):
                logger.info(" - it is full")
                return Machine(self,
                               radius=radius,
                               startangle=0.0,
                               endangle=0.0)

            logger.info("check for quarter machine")
            radius = width
            self.set_center([mm[0], mm[2]])
            logger.debug("-- center = %s,  radius min/max = %s/%s",
                         self.center, self.min_radius, self.max_radius)
            if self.check_hull(radius, mm[0], mm[2], self.rtol, atol):
                logger.info(" - it is a quarter")
                return Machine(self,
                               radius=radius,
                               startangle=0.0,
                               endangle=np.pi/2)

        elif np.isclose(width, height*2, rtol=1e-2, atol=1e-2):
            radius = width/2
            self.set_center([mm[1]-height, mm[2]])
            logger.info("check for half machine")
            if self.check_hull(radius, None, mm[2], self.rtol, atol):
                logger.info(" - it is a half")
                return Machine(self,
                               radius=radius,
                               startangle=0.0,
                               endangle=np.pi)

            self.set_center([mm[1]-height, mm[3]])
            if self.check_hull(radius, None, mm[3], self.rtol, atol):
                logger.info(" - it is a half")
                return Machine(self,
                               radius=radius,
                               startangle=np.pi,
                               endangle=0.0)

        elif np.isclose(width*2, height, rtol=1e-2, atol=1e-2):
            radius = width
            logger.info("check for half machine")
            self.set_center([mm[1], mm[3]-width])
            if self.check_hull(radius, mm[1], None, self.rtol, atol):
                logger.info(" - it is a half")
                return Machine(self, radius, np.pi/2.0, -np.pi/2.0)

            self.set_center([mm[0], mm[3]-width])
            if self.check_hull(radius, mm[0], None, self.rtol, atol):
                logger.info(" - it is a half")
                return Machine(self,
                               radius=radius,
                               startangle=-np.pi/2.0,
                               endangle=np.pi/2.0)

        machine = self.get_machine_part(mm)
        if machine:
            logger.info(" - it is 1/%s of a machine", machine.part)
            return machine

        logger.info("The shape of the Machine is unexpected")
        self.center = [0.0, 0.0]
        return Machine(self)

    def get_same_center(self, center_lst, center, rtol, atol):
        for c in center_lst:
            if points_are_close(c[1], center, rtol, atol):
                return c
        return None

    def get_machine_part(self, mm):
        logger.debug("*** Begin of get_machine_part() ***")

        h_points = [h for h in convex_hull(self.virtual_nodes())]
        center = self.get_center(h_points, mm)
        logger.debug(" - Center is %s", center)

        min_radius = 99999
        max_radius = 0
        startangle = 999.0
        endangle = -999.0

        for h in h_points:
            if not points_are_close(center, h, rtol=1e-02, atol=1e-03):
                angle = alpha_line(center, [round(h[0], 4), round(h[1], 4)])
                if angle < 0.0:
                    logger.debug(" - strange point %s", h)
                    logger.debug(" - hull-node %s ==> angle = %s", h, angle)
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

        logger.debug(" - startangle: %f12", startangle)
        logger.debug(" - endangle  : %f12", endangle)
        x = part_of_circle(startangle, endangle)
        logger.debug(" - slice     : 1/%d", x)

        logger.debug(" - center left : %f12", center_left)
        logger.debug(" - center right: %f12", center_right)
        logger.debug(" - center down : %f12", center_down)
        logger.debug(" - center up   : %f12", center_up)

        if x > 0:
            # it looks pretty good
            logger.debug(" - slice is 1/%d: EXIT", x)
            self.set_center([round(center[0], 8), round(center[1], 8)])
            M = Machine(self,
                        radius=max_radius,
                        startangle=startangle,
                        endangle=endangle)
            logger.debug("*** End of get_machine_part(): ok ***")
            return M

        if less_equal(center_left, 0.0) or less_equal(center_right, 0.0) or \
           less_equal(center_up, 0.0) or less_equal(center_down, 0.0):
            if not np.isclose(center_left, center_right):
                x = part_of_circle(startangle, endangle)

                if x > 2:
                    logger.debug(" - slice is 1/%d: EXIT", x)
                    self.set_center([round(center[0], 8), round(center[1], 8)])
                    return Machine(self,
                                   radius=max_radius,
                                   startangle=startangle,
                                   endangle=endangle)

        if min_radius >= max_radius*0.9 and min_r >= max_r*0.9:
            # Mit 10 % Abweichungen gehen wir noch von einem ganzen Motor aus.
            self.set_center([round(center[0], 8), round(center[1], 8)])
            return Machine(self,
                           radius=max(max_radius, max_r),
                           startangle=0.0,
                           endangle=0.0)

        if np.isclose(center_down, 0.0):
            min_r = min(center_left, center_right, center_up)
            if min_r >= max_r*0.9:
                # Vermutlich ein halber Motor
                self.set_center([round(center[0], 8), round(center[1], 8)])
                return Machine(self,
                               radius=max(max_radius, max_r),
                               startangle=0.0,
                               endangle=np.pi)

            if np.isclose(center_left, 0.0):
                min_r = min(center_right, center_up)
                if min_r >= max_r*0.9:
                    # Vermutlich ein viertel Motor
                    self.set_center([round(center[0], 8), round(center[1], 8)])
                    return Machine(self,
                                   radius=max(max_radius, max_r),
                                   startangle=0.0,
                                   endangle=np.pi/2)

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
                        m_min = min(m_min, m)

                y = line_n([p[0]-center[0], p[1]], m_min)
                center = (center[0], y)
                angle = alpha_line(center, p)

        self.set_center([round(center[0], 8), round(center[1], 8)])
        return Machine(self,
                       radius=max_radius,
                       startangle=angle,
                       endangle=np.pi-angle)

    def get_center_xaxis(self, points):
        y_list = [(round(p[1], 3), p) for p in points]
        y_list.sort()
        center_y = y_list[0][0]
        eq_list = [y for y, p in y_list[1:] if np.isclose(center_y, y)]
        if len(eq_list) < 2:  # min 3 Points
            return (None, None)
        start_line = Line(Element(start=(-1e+8, center_y), end=(1e+8, center_y)))
        start_m = start_line.m()
        logger.debug("Startline(%s, %s)", start_line.p1, start_line.p2)
        neq_list = [(y, p) for y, p in y_list if not np.isclose(center_y, y)]
        if len(neq_list) < 2:
            return (None, center_y)
        center_x_list = []
        for inx in range(len(neq_list)-1):
            y0, p0 = neq_list[inx]
            for y, p in neq_list[inx+1:]:
                end_line = Line(Element(start=p0, end=p))
                end_m = end_line.m()
                logger.debug("Endline(%s, %s)", end_line.p1, end_line.p2)
                if end_m is None:
                    continue
                if np.isclose(start_m, end_m):
                    continue
                p = lines_intersect_point(start_line.p1, start_m, start_line.n(start_m),
                                          end_line.p1, end_m, end_line.n(end_m))
                if p:
                    center_x_list.append(round(p[0], 2))
        center_x_list.sort()
        logger.debug("Possible center x: %s", center_x_list)
        x0 = center_x_list[0]
        c = 0
        x_list = []
        for x in center_x_list[1:]:
            if np.isclose(x0, x):
                c += 1
            else:
                x_list.append((c, x0))
                x0 = x
                c = 0
        x_list.append((c, x0))
        x_list.sort(reverse=True)
        c, center_x = x_list[0]
        if c < 3:
            center_x = None
        else:
            if (c*100 / len(neq_list)) < 10:
                center_x = None
        logger.debug("Best Center(%s, %s)", center_x, center_y)
        return (center_x, center_y)

    def get_center_hull(self, points):
        p1 = points[0]
        points.append(p1)
        p2 = points[1]
        a1 = alpha_line(p1, p2)

        logger.debug(" - p1 = %s", p1)
        logger.debug(" - p2 = %s", p2)

        lines = []
        for p in points[2:]:
            logger.debug(" - pn = %s", p)
            a2 = alpha_line(p2, p)
            if not np.isclose(a1, a2):
                d = distance(p1, p2)
                logger.debug(" - d = %s: p1/p2 = %s/%s", d, p1, p2)
                lines.append([d, a1, p1, p2])
                p1 = p2
                a1 = a2
            p2 = p

        d = distance(p1, p2)
        logger.debug(" - d = %s: p1/p2 = %s/%s", d, p1, p2)
        lines.append([d, a1, p1, p2])

        if np.isclose(lines[0][1], lines[-1][1]):
            lines[0][0] += lines[-1][0]  # distance
            lines[0][2] = lines[-1][2]   # start point
            del lines[-1]

        lines.sort(reverse=True)
        for l in lines:
            logger.debug(" - Line %s", l)

        l1 = Line(Element(start=lines[0][2], end=lines[0][3]))
        l2 = Line(Element(start=lines[1][2], end=lines[1][3]))
        center = l1.intersect_line(l2, all=True)
        logger.debug("End of get_center: %s", center)
        if center:
            return center[0]
        return None

    def is_same_center(self, center_lst, center, rtol, atol):
        for c in center_lst:
            if points_are_close(c['center'], center['center'], rtol, atol):
                radius = center['radius'][0]
                if not radius in c['radius']:
                    c['radius'].append(radius)
                c['count'] = c['count'] + 1
                return True
        return False

    def get_center_with_x(self, center_lst, x):
        c_list = [c for c in center_lst
                  if np.isclose(c['center'][0], x, rtol=self.rtol, atol=self.atol)]
        if not c_list:
            return None
        cy_list = [(c['center'][1], c) for c in c_list]
        cy_list.sort()
        y, c = cy_list[0]
        return c

    def get_center_with_y(self, center_lst, y):
        c_list = [c for c in center_lst
                  if np.isclose(c['center'][1], y, rtol=self.rtol, atol=self.atol)]
        if not c_list:
            return None
        cy_list = [(c['center'][0], c) for c in c_list]
        cy_list.sort()
        [logger.info("y=%s, c=%s", y, c) for y, c in cy_list]
        y, c = cy_list[0]
        return c

    def get_center_arcs(self, mm):
        logger.debug("begin of get_center_arcs")
        center_list = []
        x_min, x_max, y_min, y_max = mm

        def center_is_inside(c):
            tol = 0.1
            x, y = c
            if x-tol > x_min and \
               x+tol < x_max and \
               y-tol > y_min and \
               y+tol < y_max:
                return True
            return False

        circles = [e for e in self.elements() if is_Circle(e)]
        logger.debug(" -- %s Circles", len(circles))

        for e in circles:
            center = (round(e.center[0], 3), round(e.center[1], 3))
            entry = {'center': center,
                     'radius': [round(e.radius, 1)],
                     'phi': e.get_angle_of_arc(),
                     'dist': e.length(),
                     'inside': center_is_inside(center),
                     'count': 1}
            center_list.append(entry)

        arcs = [e for e in self.elements() if is_Arc(e)]
        logger.debug(" -- %s Arcs", len(arcs))

        for e in arcs:
            center = (round(e.center[0], 3), round(e.center[1], 3))
            entry = {'center': center,
                     'radius': [round(e.radius, 1)],
                     'phi': e.get_angle_of_arc(),
                     'dist': e.length(),
                     'inside': center_is_inside(center),
                     'count': 1}
            if not self.is_same_center(center_list, entry, self.rtol, self.atol):
                center_list.append(entry)

        center = None
        arc_list = [[c['count'], len(c['radius']), c['phi'], n, c]
                    for n, c in enumerate(center_list)]
        arc_list.sort(reverse=True)

        logger.debug("x min/max = %s/%s", x_min, x_max)
        logger.debug("y min/max = %s/%s", y_min, y_max)

        [logger.debug("Arc %s", arc) for arc in arc_list]
        if not arc_list:
            logger.debug("end of get_center_arcs: no arcs")
            return None

        cnt, cr1, p, n, c1 = arc_list[0]
        logger.debug("First Entry: %s", c1)
        center = c1['center']
        if len(arc_list) > 1:
            cnt, cr2, p, n, c2 = arc_list[1]
            logger.debug("Second Entry: %s", c2)
            if not cr1 > cr2:
                center = None

        if center:
            logger.debug("end of get_center_arcs: -> %s", center)
            return center

        c_entry = self.get_center_with_x(center_list, x_min)
        if c_entry:
            center = c_entry['center']
            if center[1] < y_min:
                logger.debug("end of get_center_arcs: x -> %s", center)
                return center
        c_entry = self.get_center_with_y(center_list, y_min)
        if c_entry:
            center = c_entry['center']
            if center[0] < x_min:
                logger.debug("end of get_center_arcs: y -> %s", center)
                return center

        logger.debug("end of get_center_arcs: no center found")
        return None

    def get_center_dim(self, mm):
        return (round(mm[0], 4), round(mm[2], 4))

    def get_center(self, points, mm):
        logger.debug("Begin of get_center(%s points)", len(points))
        if len(points) < 3:
            return None

        center = None
        # Zuerst suchen wir anhand der Circle- und Arc-Segmente nach einem
        # möglichen Center-Punkt.
        center_arcs = self.get_center_arcs(mm)

        if center_arcs:
            center = center_arcs
        else:
            # Wir finden keine Arc-Objekte, welche uns einen Hinweis auf den
            # Center geben können. Wir versuchen in der Verzweiflung mit
            # center_hull oder x(min) und y(min)
            center_hull = self.get_center_hull(points)
            if center_hull:
                center = center_hull
            else:
                center = self.get_center_dim(mm)

        center_xaxis = self.get_center_xaxis(points)
        y = center_xaxis[1]
        if y is not None:
            if np.isclose(y, center[1], atol=0.3):
                center = (center[0], y)
        logger.debug("End of get_center ==> %s", center)
        return center

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
                    logger.warning("BAD: Point %s", p)
                    self.airgaps.append(Point(p))
                    self.airgaps.append(e)
                    ok = False
                else:
                    borders += 1
        return (ok, borders)

    def check_airgap(self, startangle, endangle):
        logger.debug("begin check_airgap")
        area_id_list = [a.id for a in self.list_of_areas()]

        def delete_id(id):
            try:
                i = area_id_list.index(id)
            except ValueError:
                return
            area_id_list[i] = 0
        #   ---
        def append_area(alist, a):
            for my_a in alist:
                if my_a.is_in_touch_with_area(self, a):
                    alist.append(a)
                    delete_id(a.id)
                    return True
            return False
        #   ---
        for area in self.area_list:
            if not area.id in area_id_list:
                continue

            for a in self.area_list:
                if area.id == a.id:
                    continue
                if not a.id in area_id_list:
                    continue

                if area.the_area_is_inside_area(a):
                    delete_id(a.id)

        # collect remaining areas
        area_list = [a for a in self.area_list if a.id in area_id_list]
        group_list = {}
        for area in area_list:
            if not area.id in area_id_list:
                continue
            group_list[area.id] = [area]
            delete_id(area.id)
            for a in self.area_list:
                if area.id == a.id:
                    continue
                if not a.id in area_id_list:
                    continue
                if append_area(group_list[area.id], a):
                    continue

        area_list = [a for a in self.area_list if a.id in area_id_list]
        for area in area_list:
            group_list[area.id] = [area]

        area_id_list = [int(x) for x in group_list.keys()]

        logger.debug("end check_airgap: return %s", len(area_id_list) > 1)
        return len(area_id_list) > 1  # bad

    def is_border_line(self, center, startangle, endangle, e, rtol=1e-3, atol=1e-3):
        if isinstance(e, Line):
            if np.isclose(startangle, endangle):
                return False  # full

            if points_are_close(center, e.p1):
                angle_p2 = alpha_line(center, e.p2)
                if np.isclose(startangle, angle_p2, rtol=rtol, atol=atol):
                    return True
                return np.isclose(endangle, angle_p2, rtol=rtol, atol=atol)

            if points_are_close(center, e.p2):
                angle_p1 = alpha_line(center, e.p1)
                if np.isclose(startangle, angle_p1, rtol=rtol, atol=atol):
                    return True
                return np.isclose(endangle, angle_p1, rtol=rtol, atol=atol)

            angle_p1 = alpha_line(center, e.p1)
            if np.isclose(startangle, angle_p1, rtol=rtol, atol=atol):
                angle_p2 = alpha_line(center, e.p2)
                return np.isclose(startangle, angle_p2, rtol=rtol, atol=atol)

            if np.isclose(endangle, angle_p1, rtol=rtol, atol=atol):
                angle_p2 = alpha_line(center, e.p2)
                return np.isclose(endangle, angle_p2, rtol=rtol, atol=atol)
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

    def detect_airgaps(self, center, startangle, endangle, atol=0.1, with_end=False):
        """ Die Funktion sucht Luftspalt-Kandidaten und liefert eine Liste
            von Möglichkeiten mit jeweils einem minimalen und einem maximalen
            Radius als Begrenzung des Luftspalts.
        """
        gaplist = []
        for e in self.elements(Shape):
            if not self.is_border_line(center, startangle, endangle, e, atol=atol):
                gaplist += [e.minmax_from_center(center)]
            else:
                min_r, max_r = e.minmax_from_center(center)
                if with_end and np.isclose(min_r, 0.0):
                    gaplist += [(0.0, 0.0)]

        gaplist.sort()

        min_radius = gaplist[0][0] + 1
        max_radius = gaplist[-1][1] - 1
        cur_radius = gaplist[0][1]
        airgaps = []

        if with_end:
            if self.is_inner:
                min_radius = gaplist[0][0] - 1
            else:
                max_radius = gaplist[-1][1] + 1

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
                             lower_radius,
                             radius,
                             upper_radius,
                             angle_tot):
        lower_circles = self.get_circles(center, lower_radius)
        angle_sum = self.alpha_of_circles(lower_circles, center)
        if angle_sum / angle_tot < 0.5:
            return False

        upper_circles = self.get_circles(center, upper_radius)
        angle_sum = self.alpha_of_circles(upper_circles, center)
        if angle_sum / angle_tot < 0.5:
            return False

        self.remove_edges(self.get_circles(center, radius))
        return True

    def set_areas_inside_for_all_areas(self):
        logger.debug("begin set_areas_inside_for_all_areas")
        for area in self.list_of_areas():
            areas_inside = [a for a in self.area_list
                            if area.is_inside(a, self)]
            if not areas_inside:
                continue

            areas_notouch = {a.identifier(): a for a in areas_inside
                             if not area.has_connection(self, a, ndec)}
            area.areas_inside = areas_notouch

        for area in self.list_of_areas():
            logger.debug(" Inside %s is:", area.identifier())
            nested = [id for id in area.list_of_nested_areas_inside()]
            for id in nested:
                logger.debug(" -- remove %s inside %s", id, area.identifier())
                if area.areas_inside.get(id, None):
                    del area.areas_inside[id]
                else:
                    logger.warning("   %s already removed ?!", id)

        # for area in self.list_of_areas():
        #    logger.info("Inside %s is:", area.identifier())
        #    for id in area.areas_inside:
        #        logger.info(" ++ %s", id)
        logger.debug("end set_areas_inside_for_all_areas")

    def set_groups_inside_for_all_areas(self):
        logger.debug("begin set_groups_inside_for_all_areas")

        groups_inside = {}
        groups = {}
        areas_outside = []
        for area in self.list_of_areas():
            grouplist = [a for a in self.areagroup_list
                         if area.is_inside(a, self)]
            if not grouplist:
                continue
            groups_inside = {g.id: g for g in grouplist}
            area.areas_inside = groups_inside
            areas_outside.append(area)
            for g in grouplist:
                g.group_is_inside = True
                alist = groups.get(g.id, [])
                alist.append(area)
                groups[g.id] = alist

        outside_id = [a.id for a in areas_outside]
        logger.debug("Areas outside: %s", outside_id)

        for id in groups.keys():
            if len(groups[id]) > 1:
                logger.warning("Attention: nested groups of areas")
                self.journal.put("warning", "nested groups of areas")
                areas = groups[id]
                main_area = areas[0]
                main_size = main_area.area_size()
                for a in areas[1:]:
                    sz = a.area_size()
                    if sz < main_size:
                        main_area = a
                        main_size = sz
                assert(main_area is not None)
                main_area.is_child = True
                for area in areas:
                    if area.id != main_area.id:
                        del area.areas_inside[id]

        logger.debug("end set_areas_inside_for_all_areas")
        return areas_outside

    def get_minmax_magnet(self):
        logger.debug("get_minmax_magnet")
        maglist = [a for a in self.list_of_areas() if a.is_magnet()]
        min_r = 9999
        max_r = -9999
        for m in maglist:
            min_r = min(min_r, m.min_dist)
            max_r = max(max_r, m.max_dist)
        return min_r, max_r

    def create_auxiliary_lines(self, rightangle, leftangle):
        logger.debug("begin of create_auxiliary_lines")
        timer = Timer(start_it=True)
        done = False

        if True:  # new style
            logger.debug("-> start create_auxiliary_lines")
            area_list = self.list_of_areas()
            builder = AreaBuilder(geom=self)
            builder.create_area_groups(area_list)
            self.areagroup_list = builder.area_list
            area_list = self.set_groups_inside_for_all_areas()
            for area in area_list:
                if self.create_aux_lines(area, rightangle, leftangle):
                    done = True
            main_groups = [g for g in self.areagroup_list
                           if not g.group_is_inside]
            if len(main_groups) > 1:
                if self.create_outside_aux_lines(main_groups):
                    done = True
        else:
            logger.debug("-> start create_auxiliary_lines")
            self.set_areas_inside_for_all_areas()
            done = False
            for area in self.list_of_areas():
                if self.create_aux_lines(area, rightangle, leftangle):
                    done = True

        t = timer.stop("-- auxiliary lines in %0.4f seconds --")
        self.journal.put('time_auxiliary_lines', t)
        logger.debug("end of create_auxiliary_lines")
        return done

    def create_aux_lines(self, area, rightangle, leftangle):
        logger.debug("begin of create_aux_lines(%s)", area.get_id())

        areas_inside = area.areas_inside.values()
        if not areas_inside:
            logger.debug("end of create_aux_lines() for %s (no areas inside)",
                         area.get_id())
            return False

        logger.debug("areas found inside area(%s)", area.get_id())

        aux_color = 'red'
        aux_linestyle = 'dotted'
        if area.is_child:
            logger.debug("Area %s is a child of another nested area", area.id)
            rightangle = None
            leftangle = None

        areas_border = {a.get_id(): a for a in areas_inside
                        if area.has_connection(self, a, ndec)}
        areas_notouch = {a.get_id(): a for a in areas_inside
                         if not area.has_connection(self, a, ndec)}

        for id, a in areas_notouch.items():
            for id2, a2 in areas_border.items():
                if a.has_connection(self, a2, ndec):
                    areas_border[id2] = a2
                    break

        for id, a in areas_border.items():
            if id in areas_notouch:
                del areas_notouch[id]

        logger.debug("--- build notouch list ---")

        notouch_list = []
        notouch_area_list = [a for id, a in areas_notouch.items()]

        for id, a in areas_notouch.items():
            logger.debug(" --> areas_notouch: %s", id)
            if not notouch_list:
                logger.debug("   * append %s (%s)", id, a.get_id())
                notouch_list.append({id: a})
            else:
                touched = False
                for l in notouch_list:
                    for id2, a2 in l.items():
                        if a.has_connection(self, a2, ndec):
                            logger.debug("   . %s and %s are connected",
                                         a.get_id(),
                                         a2.get_id())
                            l[id] = a
                            touched = True
                            break
                    if touched:
                        break
                if not touched:
                    logger.debug("   + append %s (%s)", id, a.get_id())
                    notouch_list.append({id: a})

        logger.debug("--- search lowest gaps ---")
        created_lines = []

        # create a sorted list of lowest gaps over all areas
        # inside this area (parameter area).
        notouch_dist_list = []

        for lst in notouch_list:
            my_keys = list(lst.keys())
            my_areas = lst.values()

            low_gaps_to_parent = []
            low_gaps_to_childs = []

            for a in my_areas:
                logger.debug("==> get lowest gap from %s to %s",
                             area.get_id(),
                             a.get_id())

                low_gaps_to_parent += \
                    area.get_lowest_gap_list(a,
                                             self.center,
                                             self.max_radius,
                                             rightangle,
                                             leftangle)
                for a1 in notouch_area_list:
                    if not a1.get_id() in my_keys:
                        logger.debug("    get lowest gap from %s to %s",
                                     a.get_id(),
                                     a1.get_id())

                        low_gaps_to_childs += \
                            a1.get_lowest_gap_list(a,
                                                   self.center,
                                                   self.max_radius,
                                                   rightangle,
                                                   leftangle)

            assert(len(low_gaps_to_parent) > 0)
            low_gaps_to_parent.sort()
            low_gaps_to_childs.sort()

            aux_lines = {}
            gap_to_parent = low_gaps_to_parent[0]
            dist_to_parent = gap_to_parent[0]

            if low_gaps_to_childs:
                gap_to_child = low_gaps_to_childs[0]
                dist_to_child = gap_to_child[0]

                if dist_to_child < dist_to_parent:
                    logger.debug('gap child area is lower')
                    points = gap_to_child[1]
                    pattern = gap_to_child[2]
                    line = Line(Element(start=points[0],
                                        end=points[1]),
                                color=aux_color,
                                linestyle=aux_linestyle,
                                attr='auxline')
                    aux_lines['int'] = {'distance': dist_to_child,
                                        'p1': points[0],
                                        'p2': points[1],
                                        'line': line,
                                        'pattern': pattern,
                                        'connect': 'child',
                                        'id': gap_to_child[3]}

            my_notouch = [a for i, a in areas_notouch.items()
                          if i not in my_keys]

            points = gap_to_parent[1]
            pattern = gap_to_parent[2]
            id = gap_to_parent[3]
            connect = 'parent'
            logger.debug("   ==> try line from %s to %s",
                         points[0],
                         points[1])
            logger.debug("   Pattern is %s", pattern)
            line = Line(Element(start=points[0],
                                end=points[1]),
                        color=aux_color,
                        linestyle=aux_linestyle,
                        attr='auxline')

            intersection = False
            inner_gap_list = []
            for no_a in my_notouch:
                if no_a.intersect_area(line):
                    intersection = True
                    logger.debug("   --> intersection with %s",
                                 no_a.get_id())
                    inner_gap_list += \
                        no_a.get_lowest_gap_list(a,
                                                 self.center,
                                                 self.max_radius,
                                                 rightangle,
                                                 leftangle)

            dist = dist_to_parent
            if intersection:
                # intersection with other child area
                inner_gap_list.sort()
                dist = inner_gap_list[0][0]
                points = inner_gap_list[0][1]
                pattern = inner_gap_list[0][2]
                id = inner_gap_list[0][3]
                connect = 'child'
                logger.debug("  ==> intersection with other child area")
                line = Line(Element(start=points[0],
                                    end=points[1]),
                            color=aux_color,
                            linestyle=aux_linestyle,
                            attr='auxline')

            logger.debug("   +++ auxiliary line from %s to %s",
                         points[0], points[1])

            aux_lines['ext'] = {'distance': dist,
                                'p1': points[0],
                                'p2': points[1],
                                'line': line,
                                'pattern': pattern,
                                'connect': connect,
                                'id': id}

            if aux_lines:
                aux_lines['keys'] = my_keys
                x = len(notouch_dist_list)
                notouch_dist_list.append([dist_to_parent, x, aux_lines])

        # sort over distances to parent
        done = False
        area_to_parent_list = []
        notouch_dist_list.sort(reverse=True)

        for d, x, aux in notouch_dist_list:
            aux_line = aux.get('int', None)
            if aux_line:
                if aux_line['pattern'] in created_lines:
                    aux_line = None
            if not aux_line:
                aux_line = aux.get('ext', None)
            if aux_line:
                line = aux_line['line']
                n1 = self.find_the_node(aux_line['p1'])
                n2 = self.find_the_node(aux_line['p2'])
                logger.debug("Line: n1=%s,  n2=%s", n1, n2)

                pts = self.split_and_get_intersect_points(line)
                if len(pts) != 2:
                    logger.error("ERROR in create_aux_lines()")
                    self.journal.put("warning", "Error while creating auxiliary lines")
                    logger.debug("Points: %s", pts)
                    logger.debug("Line: %s", line)

                n1 = self.find_the_node(line.node1(ndec))
                n2 = self.find_the_node(line.node2(ndec))
                logger.debug("Line: n1=%s,  n2=%s", n1, n2)
                if n1 and n2:
                    logger.debug("Create Line %s", aux_line['pattern'])
                    self.add_element(line,
                                     rtol=self.rtol,
                                     atol=self.atol)
                    logger.debug("=== Create auxiliary line: %s", aux_line)
                    if aux_line['connect'] == 'child':
                        for a in areas_inside:
                            if a.get_id() == aux_line['id']:
                                area_to_parent_list.append(a)
                                break
                    done = True
                    created_lines.append(aux_line['pattern'])

        if not area.is_iron():
            logger.debug("end create_aux_lines() for %s (no iron)",
                         area.get_id())
            return done
        return done
        #   -----------------
        def id_of_inside_area(p):
            for a in areas_inside:
                if a.is_point_inside(p):
                    return a.id
            return 0

        #   -------------------------
        def connection_thru_main_area(pts):
            if len(pts) < 3:
                return True  # ok
            if len(areas_inside) < 2:
                return True  # ok

            id = id_of_inside_area(pts[0])
            if id == 0:  # strange
                return False  # bad

            next_id = id_of_inside_area(pts[1])
            if next_id == 0:
                return True  # ok

            id = next_id
            next_id = id_of_inside_area(pts[2])
            if id == next_id:  # thru inside-area
                return True  # ok

            return False

        #   ----------
        def takeSecond(elem):
            return elem[1]
        #   ----------

        # Arcs as additional auxiliary connections thru iron
        logger.debug("Additional Lines in %s", area.identifier())
        for a in area_to_parent_list:
            logger.debug("Area %s in Parent List", a.identifier())

            mid_dist = (a.min_dist + a.max_dist) / 2
            mid_angle = (a.min_angle + a.max_angle) / 2
            arc = Arc(Element(center=self.center,
                              radius=mid_dist,
                              start_angle=(rightangle)*180/np.pi,
                              end_angle=(mid_angle)*180/np.pi))
            pts = self.split_and_get_intersect_points(arc, aktion=False)
            pts.sort(key=takeSecond, reverse=True)
            if not connection_thru_main_area(pts):
                logger.debug("connection in nested areas")
                continue

            if len(pts) > 2:
                pts = pts[0:2]

            if len(pts) == 2:
                pts.sort(key=takeSecond)
                start_angle = alpha_line(self.center, pts[0])
                end_angle = alpha_line(self.center, pts[1])
                arc = Arc(Element(center=self.center,
                                  radius=mid_dist,
                                  start_angle=start_angle*180/np.pi,
                                  end_angle=end_angle*180/np.pi),
                          color='black',  #aux_color,
                          linestyle=aux_linestyle)
                arc.set_attribute('iron_sep')
                self.split_and_get_intersect_points(arc)
                n = self.find_nodes(pts[0], pts[1])
                self.add_or_join_edge(n[0], n[1], arc,
                                      rtol=self.rtol,
                                      atol=self.atol)

        logger.debug("end create_aux_lines() for %s",
                     area.get_id())
        return done

    def create_outside_aux_lines(self, grouplist):
        done = False
        aux_color = 'blue'
        aux_linestyle = 'dotted'

        i = 0
        gaps = []
        for group in grouplist:
            i += 1
            for g in grouplist[i:]:
                gaps += group.get_lowest_gap_list(g,
                                                  self.center, self.max_radius,
                                                  None, None)
        gaps.sort()

        l = len(grouplist) -1
        for d, points, token, id in gaps[:l]:
            logger.info("Token %s", token)
            line = Line(Element(start=points[0],
                                end=points[1]),
                                color=aux_color,
                                linestyle=aux_linestyle,
                                attr='outside_auxline')
            n1 = self.find_the_node(line.node1(ndec))
            n2 = self.find_the_node(line.node2(ndec))
            if n1 and n2:
                self.add_element(line,
                                 rtol=self.rtol,
                                 atol=self.atol)
                done = True
        return done

    def set_rotor(self):
        self.sym_counterpart = 1
        self.sym_part = 0

    def set_stator(self):
        self.sym_counterpart = 1
        self.sym_part = 2

    def force_to_be_rotor(self):
        self.sym_type = TYPE_ROTOR

    def force_to_be_stator(self):
        self.sym_type = TYPE_STATOR

    def is_rotor(self):
        if self.sym_type != TYPE_UNDEFINED:
            return self.sym_type == TYPE_ROTOR
        if self.sym_counterpart:
            return self.sym_part < self.sym_counterpart
        return False

    def is_stator(self):
        if self.sym_type != TYPE_UNDEFINED:
            return self.sym_type == TYPE_STATOR
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

        logger.debug("tiny line from %s to %s deleted", n0, n2)
        self._remove_edge(n0, n1)
        self._remove_edge(n1, n2)
        dict12['deleted'] = True
        dict01['deleted'] = True
        line = Line(Element(start=n0, end=n2))
        self.add_element(line,
                         rtol=self.rtol,
                         atol=self.atol)
        return True

    def search_tiny_elements(self, mindist):
        logger.debug("begin of search_tiny_elements(%s)", mindist)
        if mindist == 0.0:
            return []

        edges = [edge for edge in self.g.edges(data=True)
                 if distance(edge[0], edge[1]) < mindist]

        logger.debug("end of search_tiny_elements: %s tiny elements found",
                     len(edges))
        return edges

    def delete_tiny_elements(self, mindist):
        logger.debug("begin of delete_tiny_elements(%s)", mindist)
        edges = self.search_tiny_elements(mindist)
        if not edges:
            logger.debug("-- no tiny elements found")
            return 0

        deleted = 0
        for edge in edges:
            n1 = edge[0]
            n2 = edge[1]
            el = edge[2]['object']
            logger.debug("-- %s: %s <-> %s", el.classname(), n1, n2)
            logger.debug("Edge: %s", el.classname())

            if n1 is n2:
                logger.debug("-- delete edge with equal nodes")
                self._remove_edge(n1, n2)
                deleted += 1
                continue

            nbrs_n1 = [nbr for nbr in self.g.neighbors(n1)
                       if not nodes_are_equal(nbr, n2)]
            nbrs_n2 = [nbr for nbr in self.g.neighbors(n2)
                       if not nodes_are_equal(nbr, n1)]

            if len(nbrs_n1) == 0 and len(nbrs_n2) == 0:
                # lonesome edge
                logger.debug("-- delete lonesome edge")
                self._remove_edge(n1, n2)
                deleted += 1
                continue

            if len(nbrs_n1) == 1:
                if self._delete_a_tiny_element(n2, n1, edge[2], nbrs_n1[0]):
                    deleted += 1
                    continue

            if len(nbrs_n2) == 1:
                if self._delete_a_tiny_element(n1, n2, edge[2], nbrs_n2[0]):
                    deleted += 1
                    continue

        if deleted:
            logger.debug("%s tiny elements deleted", deleted)
            self.journal.put("tiny_elements_deleted", deleted)
        logger.debug("end of delete_tiny_elements")
        return deleted

    def search_critical_elements(self, mindist):
        for n in self.g.nodes():
            nbrs = self.get_neighbors(n)
            if len(nbrs) < 3:
                continue
            critical_point = False
            critical_dist = 9999
            for nbr in nbrs:
                e = self.get_edge_element(n, nbr)
                if e.is_tiny(mindist):
                    critical_point = True
                    critical_dist = min(critical_dist, e.length())
            if critical_point:
                logger.debug("Warning: maybe critical point %s", n)
                self.journal.put("maybe_critical_points", (n, critical_dist))
                c = Circle(Element(center=n, radius=1))
                self.critical_points.append(c)

    def check_shaft_area(self, shaft):
        for a in self.list_of_areas():
            if not shaft.is_identical(a):
                if shaft.is_inside(a, self):
                    shaft.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)
                    return
                if shaft.is_touching(a):
                    if not a.is_iron():
                        shaft.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)
                        return

    def mark_connecting_edges(self, windings):
        logger.debug("begin of mark_connecting_edges")
        for a in windings:
            logger.debug("- id of winding: %s", a.identifier())
            for w in windings:
                if a.id != w.id:
                    elist = [e for e in a.list_of_equal_edges(w)]
                    logger.debug(" --> %s equal egdes", len(elist))
                    for e in elist:
                        e.init_attributes('lightblue', 'no_fsl')

    def set_subregion_parameters(self,
                                 startangle,
                                 endangle):
        for a in self.list_of_areas():
            a.set_subregion_parameters(self.is_inner,
                                       self.min_radius,
                                       self.max_radius,
                                       startangle,
                                       endangle)

    def search_subregions(self, startangle, endangle, EESM, single=False):
        if self.is_stator():
            self.search_stator_subregions(startangle,
                                          endangle,
                                          single=single)
        elif self.is_rotor():
            if EESM:
                self.search_EESM_rotor_subregions(startangle,
                                                  endangle,
                                                  single=single)
            else:
                self.search_PMSM_rotor_subregions(startangle,
                                                  endangle,
                                                  single=single)
        else:
            logger.warning("no stator or rotor assigned")
            self.search_unknown_subregions()
        self.looking_for_corners()

    def get_windings(self, type):
        windings = [a for a in self.list_of_areas()
                    if a.is_type(type)]
        for w in windings:
            inside = []
            for a in self.list_of_areas():
                if not w.is_identical(a):
                    if w.is_inside(a, self):
                        inside.append(a)
            if inside:
                w.set_type(AREA.TYPE_AIR)  # air
        return [a for a in self.list_of_areas()
                if a.is_type(type)]

    def collect_windings(self):
        logger.debug("begin of collect_windings")
        good_windings = self.get_windings(AREA.TYPE_WINDINGS)
        ugly_windings = self.get_windings(AREA.TYPE_WINDINGS_OR_AIR)

        logger.debug("-- %s good and %s ugly windings",
                     len(good_windings),
                     len(ugly_windings))

        if not ugly_windings:
            logger.debug("#1 end of collect_windings: %s windings", len(good_windings))
            return good_windings

        if not good_windings:
            logger.debug("#2 end of collect_windings: %s windings", len(ugly_windings))
            [w.set_type(AREA.TYPE_WINDINGS) for w in ugly_windings]
            return [a for a in self.list_of_areas() if a.is_type(AREA.TYPE_WINDINGS)]

        # ggod and ugly windings available
        found = True
        while found:
            found = False
            for a in ugly_windings:
                if a.is_touching_areas(good_windings):
                    a.set_type(AREA.TYPE_WINDINGS)
                    found = True

            good_windings = [a for a in self.list_of_areas()
                             if a.is_type(AREA.TYPE_WINDINGS)]
            ugly_windings = [a for a in self.list_of_areas()
                             if a.is_type(AREA.TYPE_WINDINGS_OR_AIR)]

        [w.set_type(AREA.TYPE_AIR) for w in ugly_windings]
        good_windings = [a for a in self.list_of_areas()
                         if a.is_type(AREA.TYPE_WINDINGS)]

        logger.debug("return bad and ugly windings as %s good windings", len(good_windings))
        logger.debug("end of collect_windings")
        return good_windings

    def search_stator_subregions(self,
                                 startangle,
                                 endangle,
                                 single=False):
        logger.debug("Begin of search_stator_subregions")

        if self.alfa == 0.0:
            self.alfa = np.pi * 2.0

        stator_size = self.area_size()
        for area in self.list_of_areas():
            area.mark_stator_subregions(self.is_inner,
                                        stator_size,
                                        self.is_mirrored(),
                                        self.alfa,
                                        self.center,
                                        self.min_radius,
                                        self.max_radius)

        windings = self.collect_windings()
        [a.set_type(AREA.TYPE_AIR) for a in self.list_of_areas()
         if a.is_type(AREA.TYPE_WINDINGS_OR_AIR)]
        windings_found = len(windings)
        logger.debug("%d windings found", windings_found)
        self.has_windings = windings_found > 0

        if windings_found > 1:
            windings_surface = [[w.surface, w] for w in windings]
            windings_surface.sort(reverse=True)
            max_size, max_w = windings_surface[0]
            for sz, w in windings_surface[1:]:
                logger.debug("winding size = %s", sz)
                if sz / max_size < 0.70:
                    w.set_type(AREA.TYPE_AIR)
                    if sz / max_size < 0.2:
                        windings_found -= 1
            windings = [a for a in self.list_of_areas()
                        if a.is_winding()]
            if windings_found > 2 and len(windings) == 1:
                logger.info("no windings remaining")
                # no windings
                [w.set_type(AREA.TYPE_AIR) for w in windings]
                [a.set_type(AREA.TYPE_IRON) for a in self.list_of_areas()
                 if a.is_iron()]
                windings = []
            elif len(windings) < windings_found:
                logger.info("%d windings remaining", len(windings))
            if len(windings) > 2:
                self.mark_connecting_edges(windings)

        wdg_min_angle = 99999
        wdg_max_angle = 0
        wdg_min_dist = 99999
        wdg_max_dist = 0
        wdg_close_to_startangle = False
        wdg_close_to_endangle = False

        for w in windings:
            wdg_min_angle = min(wdg_min_angle, w.min_angle)
            wdg_max_angle = max(wdg_max_angle, w.max_angle)
            wdg_min_dist = min(wdg_min_dist, w.min_dist)
            wdg_max_dist = max(wdg_max_dist, w.max_dist)
            if w.close_to_startangle:
                wdg_close_to_startangle = True
            if w.close_to_endangle:
                wdg_close_to_endangle = True

        if windings_found:
            logger.debug("wdg_min_angle: %s", wdg_min_angle)
            logger.debug("wdg_max_angle: %s", wdg_max_angle)
            gap_startangle = wdg_min_angle
            gap_endangle = self.alfa - wdg_max_angle
            self.wdg_is_mirrored = self.is_mirrored()
            if np.isclose(gap_startangle,
                          gap_endangle,
                          atol=0.01):
                self.wdg_is_mirrored = False

        # air or iron near windings and near airgap ?
        air_areas = [a for a in self.list_of_areas()
                     if a.is_type(AREA.TYPE_AIR_OR_IRON)]
        for a in air_areas:
            if a.around_windings(windings, self):
                logger.debug("Area %s", a.identifier())
                logger.debug(" - air-angle min/max = %s/%s",
                             a.min_air_angle,
                             a.max_air_angle)
                logger.debug(" - wdg-angle min/max = %s/%s",
                             wdg_min_angle,
                             wdg_max_angle)

                if a.close_to_startangle:
                    if not wdg_close_to_startangle:
                        logger.debug("#0.1 ===> close to startangle")
                        a.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)
                        continue

                if a.close_to_endangle:
                    if not wdg_close_to_endangle:
                        logger.debug("#0.2 ===> close to endangle")
                        if(a.min_angle < wdg_min_angle and
                           a.close_to_ag):
                            a.set_type(AREA.TYPE_AIR)  # air
                            continue

                        a.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)
                        continue

                if greater_equal(a.min_air_angle, wdg_min_angle):
                    logger.debug("#0.3 ===> %s >= %s <===",
                                 a.min_air_angle,
                                 wdg_min_angle)

                    if a.close_to_endangle and self.is_mirrored():
                        logger.debug("#1 ===> endangle and mirrored <===")
                        a.set_type(AREA.TYPE_AIR)  # air
                    elif less_equal(a.max_air_angle, wdg_max_angle):
                        logger.debug("#2 ===> %s <= %s <===",
                                     a.max_air_angle,
                                     wdg_max_angle)
                        a.set_type(AREA.TYPE_AIR)  # air
                    else:
                        logger.debug("#3 ===> %s > %s <===",
                                     a.max_air_angle,
                                     wdg_max_angle)
                        a.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)
                else:
                    logger.debug("#4 ===> %s < %s <===",
                                 a.min_air_angle,
                                 wdg_min_angle)
                    a.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)
            else:
                logger.debug("#5 not around windings")
                a.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)

        # yoke or shaft ?
        iron_areas = [a for a in self.list_of_areas()
                      if a.is_type(AREA.TYPE_YOKE)]
        for a in iron_areas:
            if a.around_windings(windings, self):
                if less(a.min_dist, wdg_max_dist):
                    if less_equal(a.max_dist, wdg_max_dist):
                        a.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)
                    else:
                        dist_low = wdg_max_dist - a.min_dist
                        dist_up = a.max_dist - wdg_max_dist
                        if dist_low > dist_up:
                            a.set_type(AREA.TYPE_TOOTH)  # iron shaft (Zahn)

        shaft_areas = [a for a in self.list_of_areas()
                       if a.is_type(AREA.TYPE_SHAFT)]
        if shaft_areas:
            if len(shaft_areas) > 1:
                logger.debug("More than one shaft in stator ?!?")
                return
            self.check_shaft_area(shaft_areas[0])
        logger.debug("End of search_stator_subregions")

    def search_EESM_rotor_subregions(self,
                                     startangle,
                                     endangle,
                                     single=False):
        logger.debug("Begin of search_EESM_rotor_subregions")

        if self.alfa == 0.0:
            self.alfa = np.pi * 2.0

        types = {}
        for area in self.list_of_areas():
            t = area.mark_EESM_rotor_subregions(self.is_inner,
                                                self.is_mirrored(),
                                                self.alfa,
                                                self.center,
                                                self.min_radius,
                                                self.max_radius,
                                                startangle,
                                                endangle)
            if t in types:
                types[t] += 1
            else:
                types[t] = 1

        windings = [a for a in self.list_of_areas()
                    if a.is_type(AREA.TYPE_WINDINGS_OR_IRON)]

        if self.is_mirrored():
            [a.set_type(AREA.TYPE_IRON) for a in windings
             if a.close_to_endangle]
            wlist = [(a.max_angle, a) for a in self.list_of_areas()
                     if a.is_type(AREA.TYPE_WINDINGS_OR_IRON)]
            if wlist:
                wlist.sort(reverse=True)
                a, w = wlist[0]
                w.set_type(AREA.TYPE_FD_WINDINGS)
        else:
            midangle = middle_angle(startangle, endangle)
            [a.set_type(AREA.TYPE_IRON) for a in windings
             if a.max_angle > midangle and a.min_angle < midangle]
            windings = [a for a in self.list_of_areas()
                        if a.is_type(AREA.TYPE_WINDINGS_OR_IRON)]
            if len(windings) > 1:
                wlist = []
                for w in windings:
                    if w.max_angle < midangle:
                        angle = alpha_angle(w.max_angle, midangle)
                    else:
                        angle = alpha_angle(midangle, w.min_angle)
                    wlist.append((angle, w))
                wlist.sort()
                a1, w1 = wlist[0]
                a2, w2 = wlist[1]
                if np.isclose(a1, a2):
                    w1.set_type(AREA.TYPE_FD_WINDINGS)
                    w2.set_type(AREA.TYPE_FD_WINDINGS)

        # all remaining areas are in iron
        [a.set_type(AREA.TYPE_IRON) for a in self.list_of_areas()
         if a.is_type(AREA.TYPE_WINDINGS_OR_IRON)]

        logger.debug("End of search_EESM_rotor_subregions")

    def search_PMSM_rotor_subregions(self,
                                     startangle,
                                     endangle,
                                     single=False):
        logger.debug("Begin of search_PMSM_rotor_subregions")

        if self.alfa == 0.0:
            self.alfa = np.pi * 2.0

        types = {}
        for area in self.list_of_areas():
            t = area.mark_PMSM_rotor_subregions(self.is_inner,
                                                self.is_mirrored(),
                                                self.alfa,
                                                self.center,
                                                self.min_radius,
                                                self.max_radius,
                                                startangle,
                                                endangle)
            if t in types:
                types[t] += 1
            else:
                types[t] = 1

        if AREA.TYPE_SHAFT in types:
            logger.debug("Shaft is available")

        if AREA.TYPE_MAGNET_RECT_NEAR_AIRGAP in types: # magnet rectangle near ag
            if AREA.TYPE_MAGNET_RECT in types:  # magnet rectangle
                mag_rectangles = [a for a in self.list_of_areas()
                                  if a.is_type(AREA.TYPE_MAGNET_RECT)]
                for a in mag_rectangles:
                    if self.is_inner:
                        dist_1 = a.max_dist + 0.5
                        dist_2 = self.max_radius + 5
                    else:
                        dist_1 = a_min_dist - 0.5
                        dist_2 = self.min_radius - 5

                    mid_angle = a.get_mid_angle(self.center)
                    p1 = point(self.center, dist_1, mid_angle)
                    p2 = point(self.center, dist_2, mid_angle)
                    line = Line(Element(start=p1, end=p2))
                    logger.debug("magnet intersect line: %s to %s", p1, p2)
                    pts = self.split_and_get_intersect_points(line,
                                                              aktion=False,
                                                              include_end=False)
                    logger.debug("magnet intersect points: %s", pts)
                    if not pts:
                        a.set_type(AREA.TYPE_MAGNET_UNDEFINED)
                mag_rectangles = [a for a in self.list_of_areas()
                                  if a.is_type(AREA.TYPE_MAGNET_RECT)]
                if mag_rectangles:  # undo
                    logger.debug("--- undo ---")
                    [a.set_type(AREA.TYPE_MAGNET_RECT) for a in self.list_of_areas()
                     if a.is_type(AREA.TYPE_MAGNET_UNDEFINED)]
                else:
                    logger.debug("--- set type 3 ---")
                    [a.set_type(AREA.TYPE_MAGNET_AIRGAP) for a in self.list_of_areas()
                     if a.is_type(AREA.TYPE_MAGNET_RECT_NEAR_AIRGAP)]
                    [a.set_type(AREA.TYPE_MAGNET_AIRGAP) for a in self.list_of_areas()
                     if a.is_type(AREA.TYPE_MAGNET_UNDEFINED)]
                    types.pop(AREA.TYPE_MAGNET_RECT)
            else:
                # set magnet
                [a.set_type(AREA.TYPE_MAGNET_AIRGAP) for a in self.list_of_areas()
                 if a.is_type(AREA.TYPE_MAGNET_RECT_NEAR_AIRGAP)]
                types[AREA.TYPE_MAGNET_RECT_NEAR_AIRGAP] = 0

        if AREA.TYPE_MAGNET_RECT in types:  # magnet rectangle
            if types[AREA.TYPE_MAGNET_RECT] > 1:
                logger.debug("%s embedded magnets in rotor",
                             types[AREA.TYPE_MAGNET_RECT])
                emb_mag_areas = [a for a in self.list_of_areas()
                                 if a.is_type(AREA.TYPE_MAGNET_RECT)]
                [a.set_surface(self.is_mirrored()) for a in emb_mag_areas]
                max_surface = 0.0
                max_phi = 0.0
                for a in emb_mag_areas:
                    if a.surface > max_surface:
                        max_phi = a.phi
                        max_surface = a.surface

                for a in emb_mag_areas:
                    if a.surface < max_surface * 0.20:  # too small
                        logger.debug(
                            "embedded magnet %s too small: convert to air",
                            a.identifier())
                        logger.debug("max surface : %s", max_surface)
                        logger.debug("area surface: %s", a.surface)
                        logger.debug("max phi     : %s", max_phi)
                        logger.debug("area phi    : %s", a.phi)
                        if not np.isclose(a.phi, max_phi):
                            a.set_type(AREA.TYPE_AIR)  # air

            # set iron
            [a.set_type(AREA.TYPE_IRON) for a in self.list_of_areas()
             if a.is_type(AREA.TYPE_MAGNET_AIRGAP)]

        iron_mag_areas = [a for a in self.list_of_areas()
                          if a.is_type(AREA.TYPE_MAGNET_OR_IRON)]
        air_mag_areas = [a for a in self.list_of_areas()
                         if a.is_type(AREA.TYPE_MAGNET_OR_AIR)]
        ag_areas = [a for a in self.list_of_areas() if a.close_to_ag]
        if len(ag_areas) == 1:
            if len(iron_mag_areas) == 1:
                [a.set_type(AREA.TYPE_MAGNET_AIRGAP) for a in iron_mag_areas]
                iron_mag_areas = []
            if len(air_mag_areas) == 1:
                [a.set_type(AREA.TYPE_MAGNET_AIRGAP) for a in air_mag_areas]
                air_mag_areas = []

        [a.set_type(AREA.TYPE_IRON) for a in iron_mag_areas]
        [a.set_type(AREA.TYPE_AIR) for a in air_mag_areas]

        if self.is_mirrored():
            mid_alfa = round(self.alfa, 3)
        else:
            mid_alfa = round(self.alfa / 2, 4)

        mag_areas = [[abs(round(a.phi, 3) - mid_alfa),
                      a.id,
                      a] for a in self.list_of_areas()
                     if a.is_type(AREA.TYPE_MAGNET_RECT)]

        shaft_areas = [a for a in self.list_of_areas()
                       if a.is_type(AREA.TYPE_SHAFT)]
        if shaft_areas:
            if len(shaft_areas) > 1:
                logger.debug("More than one shaft in rotor ?!?")
                return
            self.check_shaft_area(shaft_areas[0])

        magnets = [a for a in self.list_of_areas()
                   if a.is_magnet()]
        self.has_magnets = len(magnets) > 0
        for m in magnets:
            m.phi = m.get_magnet_orientation()
        logger.debug("%s magnets found in rotor", len(magnets))

        if not single:
            if not magnets:
                [a.set_type(AREA.TYPE_AIR) for a in self.list_of_areas()]
                self.search_stator_subregions(startangle, endangle, single=single)
            return

        logger.debug("end of search_PMSM_rotor_subregions")

    def recalculate_magnet_orientation(self):
        logger.debug("begin of recalculate_magnet_orientation")
        magnet_areas = [a for a in self.list_of_areas()
                        if a.is_type(AREA.TYPE_MAGNET_RECT)]
        self.recalculate_magnet_group(magnet_areas)

        magnet_areas = [a for a in self.list_of_areas()
                        if a.is_type(AREA.TYPE_MAGNET_AIRGAP)]
        for a in magnet_areas:
            a.phi = a.get_magnet_orientation()
        logger.debug("end of recalculate_magnet_orientation")

    def recalculate_magnet_group(self, areas):
        if not areas:
            return
        elements = []
        for a in areas:
            elements += a.copy_of_elements()
        geom = Geometry(elements, center=self.center, type=self.sym_type)
        geom.create_list_of_areas()
        mag_areas = geom.list_of_areas()

        builder = AreaBuilder(geom=geom)
        if builder.create_area_groups(mag_areas):
            return  # bad
        group_list = builder.area_list

        for group in group_list:
            if not group.is_magnet_rectangle():
                logger.debug("Warning: group is not a rectangle")
            group.set_type(AREA.TYPE_MAGNET_RECT)
            phi = group.get_magnet_orientation()
            for a in group.areas_of_group:
                a.phi = phi
                p = a.get_best_point_inside(geom)
                phi_set = False
                for area in areas:
                    if area.the_point_is_inside_area(p):
                        logger.debug("Replace phi %s by %s", area.phi, phi)
                        area.phi = phi
                        phi_set = True
                        break
                if not phi_set:
                    logger.warning("___MAGNET PHI NOT SET. AREA NOT FOUND IN GEOMETRY___")

    def search_unknown_subregions(self):
        logger.debug("begin of search_unknown_subregions")
        types = {}
        for area in self.list_of_areas():
            t = area.mark_unknown_subregions(self.is_mirrored(),
                                             self.alfa,
                                             self.center,
                                             self.min_radius,
                                             self.max_radius)
            if t in types:
                types[t] += 1
            else:
                types[t] = 1

        if types.get(AREA.TYPE_MAGNET_RECT, 0) > 1:
            logger.debug("%s embedded magnets in rotor",
                         types[AREA.TYPE_MAGNET_RECT])
            emb_mag_areas = [a for a in self.list_of_areas()
                             if a.is_type(AREA.TYPE_MAGNET_RECT)]
            [a.set_surface(self.is_mirrored()) for a in emb_mag_areas]
            max_surface = 0.0
            for a in emb_mag_areas:
                max_surface = max(max_surface, a.surface)

            for a in emb_mag_areas:
                if a.surface < max_surface * 0.20:  # too small
                    a.set_type(AREA.TYPE_AIR)  # air

        logger.debug("end of search_unknown_subregions")

    def magnets_in_the_middle(self, midangle):
        mag_areas = [a for a in self.list_of_areas()
                     if a.is_magnet()]
        logger.debug("%s magnets in geom", len(mag_areas))
        for a in mag_areas:
            if a.max_angle > midangle and a.min_angle < midangle:
                return True
        return False

    def windings_in_the_middle(self, midangle):
        wdg_areas = [a for a in self.list_of_areas()
                     if a.is_winding()]
        logger.debug("%s windings in geom", len(wdg_areas))
        for a in wdg_areas:
            if greater(a.max_angle, midangle) and \
               less(a.min_angle, midangle):
                return True
        return False

    def looking_for_corners(self):
        if self.is_inner:
            logger.debug("looking_for_corners: inner")
            start_cp = self.start_corners[-1]
            end_cp = self.end_corners[-1]
        else:
            logger.debug("looking_for_corners: outer")
            start_cp = self.start_corners[0]
            end_cp = self.end_corners[0]
        logger.debug("looking_for_corners: start=%s, end=%s",
                     start_cp, end_cp)
        for area in self.list_of_areas():
            area.mark_airgap_corners(start_cp, end_cp)
        return

    def num_areas_of_type(self, types=()):
        return len([area for area in self.list_of_areas()
                    if area.type in types])

    def area_size_of_type(self, type):
        return sum([area.surface for area in self.list_of_areas()
                    if area.is_type(type)])*1e-3

    def num_of_windings(self):
        return self.num_areas_of_type((AREA.TYPE_WINDINGS,))

    def num_of_irons(self):
        return self.num_areas_of_type((AREA.TYPE_IRON,
                                       AREA.TYPE_YOKE,
                                       AREA.TYPE_TOOTH,))

    def num_of_magnets(self):
        return self.num_areas_of_type((AREA.TYPE_MAGNET_AIRGAP,
                                       AREA.TYPE_MAGNET_RECT,))

    def area_close_to_endangle(self, type):
        return len([area for area in self.list_of_areas()
                    if area.is_type(type) and area.close_to_endangle])

    def corners_dont_match(self):
        if self.is_mirrored():
            return False
        if len(self.start_corners) != len(self.end_corners):
            logger.warning("number of corners dont match: %s != %s",
                           len(self.start_corners),
                           len(self.end_corners))
            return True
        for i in range(len(self.start_corners)):
            d1 = distance(self.center, self.start_corners[i])
            d2 = distance(self.center, self.end_corners[i])
            if not np.isclose(d1, d2, atol=0.02):
                logger.warning("distance of corners dont match: %s != %s",
                               d1, d2)
                return True
        return False

    def get_appendices(self):
        c = 0
        end_nodes = []
        for n in self.g.nodes():
            nbrs = [nbr for nbr in self.g.neighbors(n)]
            if len(nbrs) == 1:
                el = self.get_edge_element(n, nbrs[0])
                if not el:
                    logger.error("Element not available ?!")
                    continue

                if not is_Circle(el):
                    end_nodes.append((n, nbrs[0], el))
        return end_nodes

    def connect_all_nodes(self,
                          additional_nodes=[],
                          rtol=1e-04, atol=1e-04,
                          main=False):
        logger.debug("begin of connect_all_nodes")
        timer = Timer(start_it=True)
        nodes_list = self.get_nodes()
        c1 = 0
        c2 = 0
        additional_nodes_list = []
        for n in additional_nodes:
            if not n in nodes_list:
                additional_nodes_list.append(n)
                c1 += 1
            else:
                c2 += 1
        logger.debug("connect_all_nodes: %s added, %s already available", c1, c2)

        corr = self.connect_nodes(nodes_list,
                                  rtol=rtol, atol=atol)
        if additional_nodes_list:
            corr += self.connect_nodes(additional_nodes_list,
                                       omit_single_element=True,
                                       rtol=rtol, atol=atol)
        t = timer.stop("-- {} connections in %0.4f seconds --".format(corr))
        if main:
            self.journal.put_nodes_connected(corr)
            self.journal.put('time_node_connections', t)
        logger.debug("==> %s nodes connected", corr)
        logger.debug("end of connect_all_nodes")
        return corr

    def connect_nodes(self, nodes_list,
                      rtol=1e-03, atol=1e-03,
                      omit_single_element=False,
                      ignore_end=False):
        logger.debug("begin of connect_nodes(rtol=%s, atol=%s)", rtol, atol)
        logger.debug("-- %s nodes exist", len(nodes_list))

        count = 0
        for n in nodes_list:
            if self.node_connected(n,
                                   rtol=rtol, atol=atol,
                                   omit_single_element=omit_single_element,
                                   ignore_end=ignore_end):
                count += 1

        logger.debug("end of connect_nodes => %s", count)
        return count

    def connect_all_appendices(self, rtol=1e-04, atol=1e-04, main=False):
        logger.debug("begin of connect_all_appendices")
        timer = Timer(start_it=True)
        appendix_list = self.get_appendices()
        before = len(appendix_list)
        if main:
            self.journal.put_appendices(len(appendix_list))
        corr = self.connect_appendices(appendix_list, rtol=rtol, atol=atol)
        corr_total = corr

        if corr < before:
            appendix_list = self.get_appendices()
            before = len(appendix_list)
            corr = self.connect_appendices(appendix_list, rtol=rtol, atol=atol)
            corr_total += corr
        t = timer.stop("-- {} connections in %0.4f seconds --".format(corr))
        if main:
            self.journal.put_appendices_connected(corr_total)
            self.journal.put('time_app_connections', t)

        logger.debug("==> %s appendices connected", corr_total)
        logger.debug("end of connect_all_appendices")
        return corr_total

    def connect_appendices(self, appendix_list, rtol=1e-03, atol=1e-03, ignore_end=False):
        logger.debug("begin of connect_appendices(rtol=%s, atol=%s)", rtol, atol)
        logger.debug("-- %s appendices exist", len(appendix_list))
        self.fixed_appendices = []

        count = 0
        for n0, n1, el in appendix_list:
            logger.debug("Appendix Node at %s", n0)
            if n0 in self.fixed_appendices:
                logger.debug(" - Node already fixed")
            else:
                count += self.connect_appendix(n0, n1, el,
                                               rtol=rtol, atol=atol,
                                               ignore_end=ignore_end)

        self.fixed_appendices = []

        logger.debug("end of connect_appendices => %s", count)
        return count

    def connect_appendix(self, n0, n1, el, rtol=1e-03, atol=1e-03, ignore_end=False):
        logger.debug("begin of connect_appendix(%s, rtol=%s, atol=%s)", n0, rtol, atol)

        if points_are_close(n0, n1, rtol=1e-04, atol=1e-04):
            # a very tiny appendix
            d = distance(n0, n1)
            if less(d, 0.001):
                logger.debug("-- WARNING: a very tiny appendix of length %s", d)
                nbr_list = [nbr for nbr in self.g.neighbors(n1)
                            if not nodes_are_equal(nbr, n0)]

                if len(nbr_list) == 0:
                    logger.debug("end of connect_appendix: => lonesome appendix -> no action")
                    return 0

                logger.debug("   remove it")
                self._remove_edge(n0, n1)
                return 1

        if self.node_connected(n0, rtol=rtol, atol=atol, ignore_end=ignore_end):
            logger.debug("end of connect_appendix: %s CONNECTED", n0)
            return 1

        nn = self.find_other_node(n0)
        if not nn:
            logger.debug("end of connect_appendix: => No node found nearby")
            return 0

        logger.debug("Node %s is near %s", n0, nn)
        logger.debug(" -- appendix %s from %s to %s", el.classname(), n0, n1)
        try:
            logger.debug("remove edge of %s from %s to %s",
                         el.classname(), el.p1, el.p2)
            self._remove_edge(n0, n1)
        except Exception:
            f = Path(__file__)
            msg = "{} #{}: delete of {} - {} failed".format(
                f.name, lineno(),
                n0, n1)
            self.journal.put_warning(msg)
            logger.warning("WARNING: %s", msg)
            logger.debug("-- Element %s", el)

        self.add_or_join_edge(nn, n1, el,
                              rtol=rtol,
                              atol=atol)
        self.fixed_appendices.append(nn)
        logger.debug("end of connect_appendix: connected")
        return 1

    def delete_appendices(self, appendix_list):
        c = 0
        for n0, n1, e in appendix_list:
            logger.debug("Deadend Node at %s", n0)
            c += self.remove_appendix(n0, n1)
        return c

    def delete_all_appendices(self):
        logger.debug("begin of delete_all_appendices")
        appendix_list = self.get_appendices()
        app = len(appendix_list)
        if not app:
            logger.debug("end of delete_all_appendices: no appendices")
            return

        corr = self.delete_appendices(appendix_list)
        self.journal.put_appendices_deleted(corr)

        logger.debug("==> %s appendices removed", corr)
        logger.debug("end of delete_all_appendices")

    def connect_arc_or_line(self, n, el, n1, n2, rtol=1e-03, atol=1e-03, mdec=0):
        elements = el.split([n], rtol=rtol, atol=atol, mdec=mdec)
        if len(elements) != 2:
            logger.warning("Not 2 Elements")
            logger.warning("split(rtol=%s, atol=%s, mdec=%s)", rtol, atol, mdec)
            logger.warning("Node {} in Element {}".format(n, el))
            for e in elements:
                logger.info(e)
        assert(len(elements) == 2)

        logger.debug("HIT! Node %s is in %s", n, el)
        logger.debug(" => remove %s", el)
        self.fixed_appendices.append(n1)
        self.fixed_appendices.append(n2)
        self._remove_edge(n1, n2)

        for element in elements:
            logger.debug("Split: %s", element)

        for element in elements:
            n1_inside = element.is_point_inside(n1,
                                                rtol=rtol,
                                                atol=atol,
                                                mdec=mdec,
                                                include_end=True)
            n2_inside = element.is_point_inside(n2,
                                                rtol=rtol,
                                                atol=atol,
                                                mdec=mdec,
                                                include_end=True)
            if n1_inside and n2_inside:
                logger.error("FATAL: both inside %s", element)
            elif not (n1_inside or n2_inside):
                logger.error("FATAL: neither is inside %s", element)
            else:
                if n1_inside:
                    logger.debug(" <= #1 add from %s to %s", n1, n)
                    self.add_element(element,
                                     rtol=self.rtol,
                                     atol=self.atol)
                else:
                    logger.debug(" <= #2 add from %s to %s", n2, n)
                    self.add_element(element,
                                     rtol=self.rtol,
                                     atol=self.atol)

    def connect_circle(self, n, el, n1, n2, rtol=1e-03, atol=1e-03):
        elements = el.split([n], rtol=rtol, atol=atol)
        assert(len(elements) == 3)

        logger.debug("connect_circle: Node %s is in %s", n, el)
        logger.debug(" => remove %s from %s to %s", el.classname(), n1, n2)
        e = self.get_edge_element(n1, n2)
        if not e:
            logger.error("Element from %s to %s not found", n1, n2)
        else:
            logger.debug("Element to Remove: %s", e)
        self._remove_edge(n1, n2)

        for element in elements:
            self.add_element(element,
                             rtol=self.rtol,
                             atol=self.atol)

    def node_connected(self, n,
                       rtol=1e-03, atol=1e-03,
                       omit_single_element=False,
                       ignore_end=False):
        mdec = 2
        count = 0
        el_list = []
        for n1, n2, el in self.elements_and_nodes(Shape):
            if not el.is_near(n):
                #logger.debug("Node %s is NOT near this edge", n)
                continue
            elif el.is_point_inside(n,
                                    rtol=rtol,
                                    atol=atol,
                                    mdec=mdec,
                                    include_end=False,
                                    ignore_end=ignore_end):
                logger.debug("Node %s is inside of an edge", n)
                el_list.append((n1, n2, el))

        if omit_single_element and len(el_list) < 2:
            return 0

        for n1, n2, el in el_list:
            if is_Circle(el):
                self.connect_circle(n, el, n1, n2, rtol=rtol, atol=atol)
            else:
                self.connect_arc_or_line(n, el, n1, n2, rtol=rtol, atol=atol, mdec=mdec)
            count += 1
        return count

    def remove_appendix(self, n1, n2, incr_text=''):
        e_dict = self.g.get_edge_data(n1, n2)
        if not e_dict:  # already removed
            return 0

        e = e_dict.get('object', None)
        if not e:  # unexpected
            return 0

        if isinstance(e, Circle):
            if not isinstance(e, Arc):
                return 0

        logger.debug("%s remove_appendix(%s, %s)", incr_text, n1, n2)
        self._remove_edge(n1, n2)
        c = 1
        try:
            nbrs = [nbr for nbr in self.g.neighbors(n2)]
        except nx.NetworkXError:
            logger.debug("Node %s already deleted", n2)
            nbrs = []

        if len(nbrs) == 1:
            c += self.remove_appendix(n2, nbrs[0], incr_text + '.')
        return c

    def adjust_all_points(self):
        for e in self.elements():
            e.adjust_points()

    def split_and_get_intersect_points(self, el, aktion=True, include_end=True):
        logger.debug("begin of split_and_get_intersect_points")
        rtol = 1e-03
        atol = 1e-03
        points = []
        for e in self.elements(Shape):
            pts = e.intersect_shape(el,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=include_end)
            if pts:
                logger.debug("Split %s", e)
                [logger.debug("-- intersect point %s", p) for p in pts]
                pts_inside = []
                pts_real = []
                for p in pts:
                    incl_end = is_Circle(e)
                    if not e.is_point_inside(p, rtol, atol, include_end=incl_end):
                        # get the real point
                        n = self.find_the_node(p)
                        if n:
                            for pt in points:
                                if points_are_close(pt, n):
                                    n = None
                                    break
                            if n:
                                pts_real.append(n)
                        elif e.get_point_number(p) == 1:  # strange
                            pts_real.append(e.p1)
                        else:
                            pts_real.append(e.p2)
                    else:
                        pts_real.append(p)
                        pts_inside.append(p)

                if pts_inside and aktion:
                    self.remove_edge(e)
                    elements = e.split(pts_inside, rtol, atol)
                    for e in elements:
                        n = self.find_nodes(e.start(), e.end())
                        if distance(n[0], n[1]) == 0.0:
                            logger.debug(
                                "=== OMIT ELEMENT WITH SAME NODES ===")
                        else:
                            self.add_or_join_edge(n[0], n[1], e,
                                                  rtol=rtol,
                                                  atol=atol)
                points += pts_real

        logger.debug("end of split_and_get_intersect_points")
        return points

    def _line_inside_windings(self, p1, p2):
        for area in self.list_of_areas():
            if area.is_winding():
                if area.is_point_inside(p1):
                    if area.is_point_inside(p2):
                        return True
        return False

    def _line_inside_magnets(self, p1, p2):
        for area in self.list_of_areas():
            if area.is_magnet():
                if area.is_point_inside(p1):
                    if area.is_point_inside(p2):
                        return True
        return False

    def _line_inside_air(self, p1, p2):
        for area in self.list_of_areas():
            if area.is_air():
                if area.is_point_inside(p1):
                    if area.is_point_inside(p2):
                        return True
        return False

    def _line_inside_not_iron(self, p1, p2):
        for area in self.list_of_areas():
            if area.is_shaft() or area.is_air() or area.is_magnet():
                if area.is_point_inside(p1):
                    if area.is_point_inside(p2):
                        return True
        return False

    def inside_area_list(self, p):
        for area in self.list_of_areas():
            if area.is_point_inside(p):
                yield area

    def critical_touch_point(self, points):
        logger.debug("looking for critical touch-point")
        winding_touched = False
        for p in points[1:]:
            d = distance(self.center, p)
            logger.debug("-- p = %s, dist = %s", p, d)
            for a in self.inside_area_list(p):
                logger.debug("-- Area type = %s", a.legend())
                logger.debug("        min=%s,  max= %s", a.min_dist, a.max_dist)
                logger.debug("        close to start = %s", a.close_to_startangle)
                logger.debug("        close to end   = %s", a.close_to_endangle)
                if a.is_winding():
                    winding_touched = True
                else:
                    if winding_touched and greater(a.max_dist, d, atol=0.001):
                        if not (a.close_to_startangle and a.close_to_endangle):
                            logger.debug("-- return %s", p)
                            return p
        return None

    def create_lines_outside_windings(self, points):
        if not points:
            return False
        created = False

        p1 = points[0]
        for p2 in points[1:]:
            if not points_are_close(p1, p2):
                if self._line_inside_windings(p1, p2):
                    logger.debug("in winding(%s, %s)", p1, p2)
                    p1 = p2
                    continue

                line = Line(Element(start=p1, end=p2),
                            color='darkred',
                            linestyle='dotted')
                self.add_element(line,
                                 rtol=self.rtol,
                                 atol=self.atol)
                logger.debug("add line(%s)", line)
                created = True
            p1 = p2

        return created

    def create_lines_outside_magnets(self, points):
        logger.debug("begin of create_lines_outside_magnets")
        if not points:
            return False
        created = False

        p1 = points[0]
        for p2 in points[1:]:
            logger.debug("try from %s to %s", p1, p2)
            if not points_are_close(p1, p2):
                logger.debug("Line from %s to %s", p1, p2)
                if self._line_inside_not_iron(p1, p2):
                    logger.debug("- not in iron (%s, %s)", p1, p2)
                    p1 = p2
                    continue

                line = Line(Element(start=p1, end=p2),
                            color='darkred',
                            linestyle='dotted')
                line.set_attribute('iron_sep')
                self.add_element(line,
                                 rtol=self.rtol,
                                 atol=self.atol)
                logger.debug("add line(%s)", line)
                created = True
            p1 = p2
        logger.debug("end of create_lines_outside_magnets")
        return created

    def has_areas_touching_both_sides(self):
        for a in self.area_list:
            if a.is_touching_both_sides():
                return True
        return False

    def create_inner_corner_areas(self, startangle, endangle):
        builder = AreaBuilder(geom=self)
        return builder.create_inner_corner_auxiliary_areas(startangle, endangle)

    def close_outer_winding_areas(self):
        logger.debug("begin close_outer_winding_areas(%s areas)",
                     len(self.area_list))

        builder = AreaBuilder(geom=self)
        rslt = builder.close_outer_winding_areas()
        logger.debug("end close_outer_winding_areas (%s)", rslt)
        return rslt

    def repair_border_line(self, nodes):
        logger.debug("begin repair_border_line")
        for d, n, ok in nodes:
            logger.debug(" node=%s, ok=%s", n, ok)

        d1, n1, ok1 = nodes[0]
        if not ok1:  # fatal => ignore
            logger.debug("end repair_border_line: missing point %s", n1)
            return False
        d2, n2, ok2 = nodes[-1]
        if not ok2:  # fatal => ignore
            logger.debug("end repair_border_line: missing point %s", n2)
            return False

        remove_n1 = None
        for d2, n2, ok2 in nodes[1:]:
            if ok1 and not ok2:
                remove_n1 = n1
            if ok2 and remove_n1:
                try:
                    self._remove_edge(remove_n1, n2)
                    logger.debug("Remove Line %s -- %s", remove_n1, n2)
                except nx.NetworkXError:
                    logger.debug("Warning: Remove Line %s -- %s failed", remove_n1, n2)
                    logger.debug("end repair_border_line: failed")
                    return False
                remove_n1 = None
            n1 = n2
            ok1 = ok2

        d1, n1, ok1 = nodes[0]
        for d2, n2, ok2 in nodes[1:]:
            if not ok2:  # new node
                self.add_line(n1, n2)
                logger.debug("Add Line %s -- %s", n1, n2)
            elif not ok1:
                self.add_line(n1, n2)
                logger.debug("Add Line %s -- %s", n1, n2)
            n1 = n2
            ok1 = ok2
        logger.debug("end repair_border_line")
        return True

    def create_boundary_nodes(self,
                              center,
                              startangle,
                              endangle,
                              rtol=None, atol=None):
        logger.debug("begin of create_boundary_nodes")
        if not rtol:
            rtol = 1e-3
        if not atol:
            atol = 1e-3

        def check_line(nlist):
            d, n1, b = nlist[0]
            for d, n2, b in nlist[1:]:
                if not self.get_edge_element(n1, n2):
                    return False
                n1 = n2
            return True

        start_nodes = [(distance(center, n), n, True)
                       for n in self.angle_nodes(center, startangle, rtol, atol)]
        start_nodes.sort()
        d_start1, n, b = start_nodes[0]
        if not points_are_close(self.start_corners[0], n):
            logger.warning("end of create_boundary_nodes: corner missing in start boundary")
            return False
        d_start2, n, b = start_nodes[-1]
        if not points_are_close(self.start_corners[-1], n):
            logger.warning("end of create_boundary_nodes: corner missing in start boundary")
            return False
        if not check_line(start_nodes):
            logger.warning("end of create_boundary_nodes: bad start boundary")
            return False

        logger.debug("Start Nodes")
        [logger.debug(" --> %s", x) for x in start_nodes]

        end_nodes = [(distance(center, n), n, True)
                     for n in self.angle_nodes(center, endangle, rtol, atol)]
        end_nodes.sort()
        d_end1, n, b = end_nodes[0]
        if not points_are_close(self.end_corners[0], n):
            logger.warning("end of create_boundary_nodes: corner missing in end boundary")
            return False
        d_end2, n, b = end_nodes[-1]
        if not points_are_close(self.end_corners[-1], n):
            logger.warning("end of create_boundary_nodes: corner missing in end boundary")
            return False
        if not check_line(end_nodes):
            logger.warning("end of create_boundary_nodes: bad end boundary")
            return False

        logger.debug("End Nodes")
        [logger.debug(" --> %s", x) for x in end_nodes]

        logger.debug("Lower Corners: %s <> %s", d_start1, d_end1)
        if not np.isclose(d_start1, d_end1, rtol=self.rtol, atol=self.atol):
            logger.warning("end of create_boundary_nodes: corners dont match")
            return False

        logger.debug("Upper Corners: %s <> %s", d_start2, d_end2)
        if not np.isclose(d_start2, d_end2, rtol=self.rtol, atol=self.atol):
            logger.warning("end of create_boundary_nodes: corners dont match")
            return False

        if len(start_nodes) == 2 and len(end_nodes) == 2:
            logger.debug("end of create_boundary_nodes: only corners available")
            return False  # ok

        def node_distance_list(nodelist1, nodelist2):
            distlist = []
            i1 = 0
            i2 = 0
            while i1 < len(nodelist1) and i2 < len(nodelist2):
                d1, n1, b1 = nodelist1[i1]
                d2, n2, b2 = nodelist2[i2]
                if np.isclose(d1, d2, rtol=self.rtol, atol=self.atol):
                    distlist.append((d1, True, True))
                    i1 += 1
                    i2 += 1
                elif d1 > d2:
                    distlist.append((d2, False, True))
                    i2 += 1
                else:
                    distlist.append((d1, True, False))
                    i1 += 1
            if not i1 == len(nodelist1) and i2 == len(nodelist2):
                return []
            return distlist

        distance_list = node_distance_list(start_nodes, end_nodes)
        [logger.debug("distance: %s, (%s, %s)", d, b1, b2)
         for d, b1, b2 in distance_list]

        diff = len(distance_list) - len(start_nodes)
        done = False
        if not diff == 0:
            logger.debug("%s missing start nodes", diff)
            done = True
            for d, in_start, in_end in distance_list:
                if not in_start:
                    p = point(self.center, d, startangle)
                    start_nodes.append((d, p, False))
            start_nodes.sort()
            assert(len(start_nodes) == len(distance_list))
            if not self.repair_border_line(start_nodes):
                logger.debug("end of create_boundary_nodes (failed)")
                return False

        diff = len(distance_list) - len(end_nodes)
        if not diff == 0:
            logger.debug("%s missing end nodes", diff)
            done = True
            for d, in_start, in_end in distance_list:
                if not in_end:
                    p = point(self.center, d, endangle)
                    end_nodes.append((d, p, False))
            end_nodes.sort()
            assert(len(end_nodes) == len(distance_list))

            if not self.repair_border_line(end_nodes):
                logger.debug("end of create_boundary_nodes (failed)")
                return False

        logger.debug("end of create_boundary_nodes")
        return done

    def set_point_of_node(self, node, p):
        if isinstance(node, list):
            node = (node[0], node[1])
        nbrs = self.get_neighbors(node)
        for nbr in nbrs:
            e = self.get_edge_element(node, nbr)
            e.replace_point(node, p)

    def create_and_append_area(self, n1, n2):
        rslt = self.get_new_area(n1, n2, False)
        logger.debug("create_and_append_area: %s", rslt)
        if rslt.get('ok', False):
            area = rslt['area']
            a = Area(area, self.center, 0.0)
            a.set_type(AREA.TYPE_AIR)  # air
            self.area_list.append(a)
            return True
        logger.error("No area for air near airgap!!")
        return False

    def get_intersection_points(self, elements, line, n):
        points = []
        for e in elements:
            pts = e.intersect_line(line, include_end=True)
            if len(pts) == 1:
                if points_are_close(pts[0], n):
                    continue
            points += pts
        return points

    def mirror_all_areas(self, mirrorangle):
        axis_p = point(self.center, self.max_radius, mirrorangle)
        axis_m = line_m(self.center, axis_p)
        axis_n = line_n(self.center, axis_m)

        def add_element(e):
            n = self.find_nodes(e.start(), e.end())
            if not self.has_edge(n[0], n[1]):
                self.add_edge(n[0], n[1], e)

        area_list = []
        for a in self.area_list:
            area = a.mirror_area(self.center, axis_m, axis_n)
            area_list.append(area)
            for e in area.elements():
                add_element(e)
        self.area_list += area_list

    def areas_minmax_list(self, area_list):
        if not area_list:
            return []

        dist_list = []
        for n, a in enumerate(area_list):
            dist_list.append((a.min_dist, a.max_dist, n))
        dist_list.sort()

        minmax_list = []
        d1_this, d2_this, n = dist_list [0]
        for d1_next, d2_next, n in dist_list[1:]:
            if d1_next > d2_this:
                minmax_list.append((d1_this, d2_this))
                d1_this = d1_next
                d2_this = d2_next
            else:
                d2_this = max(d2_this, d2_next)
        minmax_list.append((d1_this, d2_this))
        return minmax_list

    def magnets_minmax_list(self):
        magnets = [a for a in self.list_of_areas()
                   if a.is_magnet()]
        return self.areas_minmax_list(magnets)

    def windings_minmax_list(self):
        windings = [a for a in self.list_of_areas()
                    if a.is_winding()]
        return self.areas_minmax_list(windings)

    def fd_windings_minmax_list(self):
        windings = [a for a in self.list_of_areas()
                    if a.is_field_winding()]
        return self.areas_minmax_list(windings)

    def shaft_minmax(self):
        shafts = [a for a in self.list_of_areas()
                  if a.is_shaft()]
        if not shafts:
            return (0.0, 0.0)

        dmin = shafts[0].min_dist
        dmax = shafts[0].max_dist
        for s in shafts[1:]:
            dmin = min(dmin, s.min_dist)
            dmax = max(dmax, s.max_dist)
        return (dmin, dmax)

    def check_airgap_connecting_nodes(self, geom, startangle, endangle):
        logger.info("check_airgap_connecting_nodes")

        start_node_in = self.get_start_airgap_node()
        start_node_out =  geom.get_start_airgap_node()
        angle_in = alpha_line(self.center, start_node_in)
        angle_out = alpha_line(self.center, start_node_out)
        if np.isclose(angle_in, angle_out):
            return
        logger.warning("WARNING: airgap connecting nodes do not aline")

    def print_nodes(self):
        print("=== List of Nodes ({}) ===".format(self.number_of_nodes()))
        for n in self.g.nodes():
            print(n)

    def write_nodes(self, filename):
        global nodes_filecount
        ts = time.localtime(time.time())
        nodes_filecount += 1
        name = "{}_{}_{:02}_{:02}_{:02}.{}".format(filename,
                                                   nodes_filecount,
                                                   ts.tm_hour,
                                                   ts.tm_min,
                                                   ts.tm_sec,
                                                   'txt')
        logger.info(">>> write nodes %s <<<", name)
        nodes = []
        for n in self.g.nodes():
            nodes.append(n)
        nodes.sort()

        content = []
        for n in nodes:
            content.append(u"{}".format(n))
            for nbr in self.g.neighbors(n):
                e_dict = self.g.get_edge_data(n, nbr)
                if not e_dict:
                    raise ValueError("Fatal: no edge-data found")
                e = e_dict['object']
                content.append(u" --> {}   e-nodes {} / {}"
                               .format(nbr, e.n1, e.n2))

        with io.open(name, 'w', encoding='utf-8') as f:
            f.write('\n'.join(content))
