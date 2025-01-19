# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.shape
  ~~~~~~~~~~~~~~~~~~~~~~

  shape objects are basic geometric elements which have 2 nodes

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function
import numpy as np
import logging
from .functions import less_equal, greater_equal
from .functions import distance, line_m, line_n, mirror_point
from .functions import point, points_are_close, points_on_arc
from .functions import alpha_line, alpha_angle, alpha_triangle
from .functions import normalise_angle, min_angle, max_angle, get_angle_of_arc
from .functions import lines_intersect_point, nodes_are_equal
from .functions import is_angle_inside, intersect_point, point_greater_equal
from .functions import middle_angle, middle_point_of_line, elevation_angle

logger = logging.getLogger('femagtools.geom')


class Element(object):
    """value object class"""

    def __init__(self, **kwargs):
        for k in kwargs.keys():
            setattr(self, k, kwargs[k])


class Shape(object):
    def init_attributes(self, color=None, attr=None, linestyle=None):
        if color is not None:
            self.my_color = color
        if linestyle is not None:
            self.my_linestyle = linestyle
        if attr is not None:
            if not hasattr(self, 'my_attrs'):
                self.my_attrs = []
            self.my_attrs.append(attr)

    def copy_attributes(self, s):
        if hasattr(s, 'my_color'):
            self.my_color = s.my_color
        if hasattr(s, 'my_linestyle'):
            self.my_linestyle = s.my_linestyle
        if hasattr(s, 'my_attrs'):
            self.my_attrs = []
            for a in s.my_attrs:
                self.my_attrs.append(a)

    def classname(self):
        return "Shape"

    def clone(self):
        return None

    def set_my_color(self, color):
        self.my_color = color

    def get_my_color(self):
        if hasattr(self, 'my_color'):
            return self.my_color
        return None

    def get_my_linestyle(self):
        if hasattr(self, 'my_linestyle'):
            return self.my_linestyle
        return None

    def set_attribute(self, attr):
        if attr is None:
            return
        if not hasattr(self, 'my_attrs'):
            self.my_attrs = []
        self.my_attrs.append(attr)

    def has_attribute(self, attr):
        if hasattr(self, 'my_attrs'):
            return attr in self.my_attrs
        return False

    def set_nodes(self, n1, n2):
        if self.get_point_number(n1) == 1:
            self.n1 = n1
            self.n2 = n2
        else:
            self.n2 = n1
            self.n1 = n2

    """an abstract geometry with 2 points"""

    def start(self):
        return self.p1

    def end(self):
        return self.p2

    def node1(self, ndec):
        return round(self.p1[0], ndec), round(self.p1[1], ndec)

    def node2(self, ndec):
        return round(self.p2[0], ndec), round(self.p2[1], ndec)

    def xmin(self):
        return min(self.p1[0], self.p2[0])

    def xmax(self):
        return max(self.p1[0], self.p2[0])

    def ymin(self):
        return min(self.p1[1], self.p2[1])

    def ymax(self):
        return max(self.p1[1], self.p2[1])

    def dx(self):
        return (self.p2[0]-self.p1[0])

    def dy(self):
        return (self.p2[1]-self.p1[1])

    def m(self, none_val=None, dec=0):
        return line_m(self.p1, self.p2, none_val=none_val, dec=dec)

    def n(self, m):
        return line_n(self.p1, m)

    def move(self, dist, ndec):
        self.p1 = self.p1[0] + dist[0], self.p1[1] + dist[1]
        self.p2 = self.p2[0] + dist[0], self.p2[1] + dist[1]
        self.n1 = (round(self.n1[0] + dist[0], ndec),
                   round(self.n1[1] + dist[1], ndec))
        self.n2 = (round(self.n2[0] + dist[0], ndec),
                   round(self.n2[1] + dist[1], ndec))

    def scale(self, factor, ndec):
        self.p1 = factor*self.p1[0], factor*self.p1[1]
        self.p2 = factor*self.p2[0], factor*self.p2[1]
        self.n1 = (round(factor*self.p1[0], ndec),
                   round(factor*self.p1[1], ndec))
        self.n2 = (round(factor*self.p2[0], ndec),
                   round(factor*self.p2[1], ndec))

    def transform(self, T, alpha, ndec):
        n = T.dot(np.array((self.p1[0], self.p1[1])))
        self.p1 = (n[0], n[1])
        n = T.dot(np.array((self.p2[0], self.p2[1])))
        self.p2 = (n[0], n[1])
        if not self.n1:
            return self
        n = T.dot(np.array((self.n1[0], self.n1[1])))
        self.n1 = (round(n[0], ndec),
                   round(n[1], ndec))
        n = T.dot(np.array((self.n2[0], self.n2[1])))
        self.n2 = (round(n[0], ndec),
                   round(n[1], ndec))
        return self

    def correct(self, src_alpha, dest_alpha, ndec):
        alpha_ref = min(src_alpha, dest_alpha)
        alpha_p = alpha_line((0.0, 0.0), self.n1)
        if greater_equal(alpha_p, alpha_ref):
            delta = dest_alpha - alpha_p
            T = np.array(((np.cos(delta), -np.sin(delta)),
                          (np.sin(delta), np.cos(delta))))
            n = T.dot(np.array((self.p1[0], self.p1[1])))
            self.p1 = (n[0], n[1])
            n = T.dot(np.array((self.n1[0], self.n1[1])))
            self.n1 = (round(n[0], ndec),
                       round(n[1], ndec))

        alpha_p = alpha_line((0.0, 0.0), self.n2)
        if greater_equal(alpha_p, alpha_ref):
            delta = dest_alpha - alpha_p
            T = np.array(((np.cos(delta), -np.sin(delta)),
                          (np.sin(delta), np.cos(delta))))
            n = T.dot(np.array((self.p2[0], self.p2[1])))
            self.p2 = (n[0], n[1])
            n = T.dot(np.array((self.n2[0], self.n2[1])))
            self.n2 = (round(n[0], ndec),
                       round(n[1], ndec))

        return self

    def overlapping_shapes(self, n, e, rtol=1e-03, atol=1e-03):
        return False

    def overlapping_shape(self, e, rtol=1e-03, atol=1e-03):
        # element e is already installed
        return None  # no overlap

    def intersect_shape(self, e, rtol=1e-03, atol=1e-03, include_end=False):
        if isinstance(e, Line):
            return self.intersect_line(e, rtol, atol, include_end)
        if isinstance(e, Arc):
            return self.intersect_arc(e, rtol, atol, include_end)
        if isinstance(e, Circle):
            return self.intersect_circle(e, rtol, atol, include_end)
        return []

    def get_node_number(self, n, override=False):
        if not nodes_are_equal(n, self.n1):
            if not nodes_are_equal(n, self.n2):
                if override:
                    return 0
                logger.debug("FATAL: get_node_number(): node %s missing in %s",
                             n, self)
                raise ValueError('missing node in element')

        d1 = distance(n, self.n1)
        d2 = distance(n, self.n2)
        if d1 == d2:
            logger.warning("distances of %s and %s are equal (%s / %s)",
                           self.n1, self.n2, d1, d2)
            #raise ValueError('both nodes are equal in element')
            logger.warning("Points of %s are %s and %s",
                           self.classname(), self.p1, self.p2)
            return 0

        if d1 < d2:
            return 1
        else:
            return 2

    def get_point_number(self, p):
        d_p1 = distance(p, self.p1)
        d_p2 = distance(p, self.p2)
        if d_p1 < d_p2:
            return 1
        else:
            return 2

    def get_alpha(self, n):
        return 0.0

    def minmax_angle_dist_from_center(self,
                                      min_alfa,
                                      max_alfa,
                                      center,
                                      circ):
        points = self.intersect_circle(circ, include_end=True)
        if not points:
            return None
        my_min_angle = min_alfa
        my_max_angle = max_alfa
        for p in points:
            alpha = alpha_line(center, p)
            my_min_angle = min_angle(my_min_angle, alpha)
            my_max_angle = max_angle(my_max_angle, alpha)

        return (my_min_angle, my_max_angle)

    def concatenate(self, n1, n2, el,
                    rtol=1e-05, atol=1e-05,
                    mdec=0,
                    overlapping=False):
        if isinstance(el, Line):
            return self.concatenate_line(n1, n2, el,
                                         rtol=rtol, atol=atol,
                                         mdec=mdec,
                                         overlapping=overlapping)
        if isinstance(el, Arc):
            return self.concatenate_arc(n1, n2, el,
                                        rtol=rtol, atol=atol,
                                        overlapping=overlapping)

        if isinstance(el, Circle):
            return self.concatenate_circle(n1, n2, el,
                                           rtol=rtol, atol=atol,
                                           overlapping=overlapping)
        return None

    def concatenate_line(self, n1, n2, el,
                         rtol=1e-05, atol=1e-05,
                         mdec=0,
                         overlapping=False):
        return None

    def concatenate_arc(self, n1, n2, el,
                        rtol=1e-05, atol=1e-05,
                        overlapping=False):
        return None

    def concatenate_circle(self, n1, n2, el,
                           rtol=1e-05, atol=1e-05,
                           overlapping=False):
        return None

    def points_sorted(self, rtol=1e-05, atol=1e-05):
        if point_greater_equal(self.p1, self.p2, rtol=rtol, atol=atol):
            return self.p2, self.p1
        else:
            return self.p1, self.p2

    def rotate(self, T, p):
        n = T.dot(np.array((p[0], p[1])))
        return (n[0], n[1])

    def is_tiny(self, mindist):
        return distance(self.p1, self.p2) < mindist

    def adjust_points(self):
        self.p1 = self.n1
        self.p2 = self.n2

    def replace_point(self, node, point):
        if nodes_are_equal(node, self.n1):
            self.p1 = point
        elif nodes_are_equal(node, self.n2):
             self.p2 = point
        else:
            logger.warning("Node %s not in element", node)

    def is_near(self, n):
        return False

    def mirror_shape(self, geom_center, axis_m, axis_n):
        n2 = mirror_point(self.start(), geom_center, axis_m, axis_n)
        n1 = mirror_point(self.end(), geom_center, axis_m, axis_n)

        el = None
        if isinstance(self, Line):
            el = Line(Element(start=n1, end=n2))

        elif isinstance(self, Arc):
            c = mirror_point(self.center, geom_center, axis_m, axis_n)
            alpha1 = alpha_line(c, n1)
            alpha2 = alpha_line(c, n2)
            el = Arc(Element(center=c,
                             radius=self.radius,
                             start_angle=alpha1*180/np.pi,
                             end_angle=alpha2*180/np.pi))

        elif isinstance(self, Circle):
            c = mirror_point(self.center, geom_center, axis_m, axis_n)
            el = Circle(Element(center=c,
                                radius=self.radius))

        if el:
            el.copy_attributes(self)
        return el

    def print_nodes(self):
        return " n1={}/n2={}".format(self.n1, self.n2)

    def __str__(self):
        return " {}/{}".format(self.p1, self.p2)

    def __lt__(self, s):
        return False


#############################
#       Circle (Shape)      #
#############################

class Circle(Shape):
    """a circle with center and radius"""

    def __init__(self, e, lf=1,
                 color=None, attr=None, linestyle=None,
                 xoff=0.0, yoff=0.0, rotation=0.0):
        self.init_attributes(color, attr, linestyle)
        if rotation != 0.0:
            alpha = rotation * np.pi/180
            T = np.array(((np.cos(alpha), -np.sin(alpha)),
                          (np.sin(alpha), np.cos(alpha))))
            center = self.rotate(T, e.center)
        else:
            center = e.center
        self.center = lf*center[0] + xoff, lf*center[1] + yoff
        self.radius = lf*e.radius
        self.startangle = 0.0
        self.endangle = 0.0
        self.p1 = self.center[0]-self.radius, self.center[1]
        self.p2 = self.center[0]+self.radius, self.center[1]
        self.n1 = None
        self.n2 = None

    def classname(self):
        return "Circle"

    def clone(self):
        return Circle(Element(center=self.center,
                              radius=self.radius))

    def render(self, renderer, color='blue', with_nodes=False):
        tmp_color = self.get_my_color()
        if not tmp_color:
            tmp_color = color
        tmp_linestyle = self.get_my_linestyle()

        renderer.circle(self.center, self.radius,
                        color=tmp_color,
                        linestyle=tmp_linestyle)
        if with_nodes:
            renderer.point(self.center, 'ro', 'white')

    def move(self, dist, ndec):
        super(Circle, self).move(dist, ndec)
        self.center = self.center[0]+dist[0], self.center[1]+dist[1]

    def minmax(self):
        """ Die Funktion bestimmt das Minimum und Maximum auf der x- und der
            y-Achse (return [<min-x>, <max-x>, <min-y>, <max-y>])
        """
        return [self.center[0]-self.radius, self.center[0]+self.radius,
                self.center[1]-self.radius, self.center[1]+self.radius]

    def minmax_from_center(self, center):
        """ Die Funktion ermittelt den minimalen und maximalen Abstand vom Center
        """
        d = distance(center, self.center)
        if np.isclose(d, 0.0):
            return (self.radius, self.radius)

        dist_min = abs(d - self.radius)
        dist_max = d + self.radius
        return (dist_min, dist_max)

    def minmax_angle_from_center(self, center):
        d = distance(center, self.center)
        r = self.radius
        if r >= d:
            return (0.0, 0.0)

        r2 = np.sqrt(d**2 - r**2)
        circ = Circle(Element(center=center, radius=r2))
        points = self.intersect_circle(circ)

        if len(points) == 2:
            alpha_p1 = alpha_line(center, points[0])
            alpha_p2 = alpha_line(center, points[1])
            if alpha_angle(alpha_p1, alpha_p2) < np.pi:
                return (alpha_p1, alpha_p2)
            else:
                return (alpha_p2, alpha_p1)
        else:
            return (0.0, 0.0)

    def get_nodes(self, parts=8, render=False):
        """ returns a list of virtual nodes to create the convex hull
        """
        return (p for p in points_on_arc(self.center, self.radius,
                                         0.0,  # startangle
                                         0.0,  # endangle
                                         parts=parts))

    def scale(self, factor):
        super(Circle, self).scale(factor)
        self.center = factor*self.center[0], factor*self.center[1]
        self.radius = factor*self.radius

    def transform(self, T, alpha, ndec):
        super(Circle, self).transform(T, alpha, ndec)
        n = T.dot(np.array((self.center[0], self.center[1])))
        self.center = (n[0], n[1])
        return self

    def center_of_connection(self, ndec=6):
        return (self.center[0] + self.radius, self.center[1])

    def maxdist(self, r):
        r0 = np.linalg.norm(self.center)
        return max(abs(r-r0-self.radius), abs(r-r0+self.radius))

    def length(self):
        """returns length of this circle"""
        return self.radius*2*np.pi

    def overlapping_shape(self, e, rtol=1e-03, atol=1e-03):
        if not (isinstance(e, Arc) or isinstance(e, Circle)):
            # end overlapping_shape: Circle (not Arc or Circle)
            return None

        if not (points_are_close(self.center, e.center) and
                np.isclose(self.radius, e.radius)):
            # end overlapping_shape: Circle (radius/center are not equal)
            return None

        # overlapping candidates
        if not isinstance(e, Arc):  # it's a circle
            # end overlapping_shape: Circle (OVERLAPPING CIRCLES)
            return [e]  # The same circle twice

        # Replace CIRCLE by two Arcs
        points = []
        points.append(e.p1)
        points.append(e.p2)
        points.append(e.p1)
        return self.create_arcs(points)

    def create_arcs(self, points):
        if not points:
            return None

        pieces = []
        p1 = points[0]
        for p2 in points[1:]:
            alpha1 = alpha_line(self.center, p1)
            alpha2 = alpha_line(self.center, p2)
            arc = Arc(Element(center=self.center,
                              radius=self.radius,
                              start_angle=alpha1*180/np.pi,
                              end_angle=alpha2*180/np.pi))
            pieces.append(arc)
            p1 = p2
        return pieces

    def intersect_line(self, line, rtol=1e-03, atol=1e-03, include_end=False):
        """ Von einem Circle-Objekt und einem Line-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben.
        """
        line_m = line.m()
        p = []
        if line_m is None:
            p = [line.p1[0], self.center[1]]
        elif np.isclose(line_m, 0.0, rtol, atol):
            p = [self.center[0], line.p1[1]]
        else:
            m = -1/line_m
            p = lines_intersect_point(line.p1, line_m, line.n(line_m),
                                      self.center, m, line_n(self.center, m))

        d = distance(self.center, p)

        if np.isclose(d, self.radius, rtol, atol):
            if line.is_point_inside(p,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=include_end):
                # Wenn der Abstand d dem Radius entspricht, handelt es sich um
                # eine Tangente und es gibt genau einen Schnittpunkt
                if include_end:
                    return [p]
                else:
                    return []
        if self.radius < d:
            # d liegt ausserhalb des Kreises -> kein Schnittpunkt
            return []

        A = np.sqrt(self.radius**2 - d**2)
        delta = alpha_line(line.p1, line.p2)

        p1 = point(p, -A, delta)
        p2 = point(p, A, delta)

        # Die Schnittpunkte p1 und p2 sind bestimmt. Nun muss noch sicher
        # gestellt werden, dass sie innerhalb des Start- und Endpunkts der
        # Linie liegen
        p1_inside = line.is_point_inside(p1,
                                         rtol=rtol,
                                         atol=atol,
                                         include_end=include_end)
        p2_inside = line.is_point_inside(p2,
                                         rtol=rtol,
                                         atol=atol,
                                         include_end=include_end)
        if p1_inside:
            if p2_inside:
                return [p1, p2]
            else:
                return[p1]
        else:
            if p2_inside:
                return[p2]
            else:
                return []

    def intersect_circle(self, circle, rtol=1e-03, atol=1e-03,
                         include_end=False):
        """ Von zwei Circle-Objekten werden die Schnittpunkte bestimmt
            und in einer Liste ausgegeben
        """
        d = distance(self.center, circle.center)
        arc = alpha_triangle(circle.radius, self.radius, d)

        if np.isnan(arc):
            if not np.isclose(d, circle.radius + self.radius,
                              rtol, atol):
                return []
            arc = 0.0
        arc_C = alpha_line(self.center, circle.center)
        p1 = point(self.center, self.radius, arc_C+arc)
        p2 = point(self.center, self.radius, arc_C-arc)
        if points_are_close(p1, p2, rtol, atol):
            # Tangente
            if include_end:
                return [p1]
            else:
                return []
        return [p1, p2]

    def intersect_arc(self, arc, rtol=1e-03, atol=1e-03, include_end=False):
        """ Von einem Circle-Objekt und einem Arc-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben
        """
        assert(isinstance(arc, Arc))
        # let Arc do the work
        return arc.intersect_circle(self, rtol, atol, include_end)

    def concatenate_arc(self, n1, n2, el,
                        rtol=1e-03, atol=1e-03,
                        overlapping=False):
        if not points_are_close(self.center, el.center):
            return None
        if not np.isclose(self.radius, el.radius):
            return None
        # it's a circle
        return Circle(Element(center=self.center, radius=self.radius))

    def concatenate_circle(self, n1, n2, el,
                           rtol=1e-03, atol=1e-03,
                           overlapping=False):
        if not points_are_close(self.center, el.center):
            return None
        if not np.isclose(self.radius, el.radius):
            return None
        # it's a circle
        return Circle(Element(center=self.center, radius=self.radius))

    def is_point_inside(self, p,
                        rtol=1e-03,
                        atol=1e-03,
                        include_end=False,
                        ignore_end=False,
                        mdec=0):
        """ returns true if p is on circle
        """
        d = distance(p, self.center)
        if not np.isclose(d, self.radius, rtol=rtol, atol=atol):
            return False
        if points_are_close(p, self.p1, rtol=rtol, atol=atol):
            return include_end
        elif points_are_close(p, self.p2, rtol=rtol, atol=atol):
            return include_end
        return True

    def split(self, points, rtol=1e-03, atol=1e-03, mdec=0):
        """ Die Funktion splittet das Circle-Objekt an den vorgegebenen Punkten
            und gibt eine Liste der neu enstandenen Elemente aus.
        """
        if len(points):
            p = points[0]
            split_arcs = []
            alpha1 = alpha_line(self.center, p)
            if len(points) == 2:
                p = points[1]
                alpha3 = alpha_line(self.center, p)
                alpha2 = middle_angle(alpha1, alpha3)
                alpha4 = middle_angle(alpha3, alpha1)
            else:
                alpha2 = normalise_angle(alpha1 + np.pi/2)
                alpha3 = normalise_angle(alpha1 - np.pi/2)
                alpha4 = None

            arc = Arc(Element(center=self.center, radius=self.radius,
                              start_angle=alpha1*180/np.pi,
                              end_angle=alpha2*180/np.pi))
            arc.copy_attributes(self)
            split_arcs.append(arc)

            arc = Arc(Element(center=self.center, radius=self.radius,
                              start_angle=alpha2*180/np.pi,
                              end_angle=alpha3*180/np.pi))
            arc.copy_attributes(self)
            split_arcs.append(arc)

            if alpha4:
                arc = Arc(Element(center=self.center, radius=self.radius,
                                  start_angle=alpha3*180/np.pi,
                                  end_angle=alpha4*180/np.pi))
                arc.copy_attributes(self)
                split_arcs.append(arc)
            else:
                alpha4 = alpha3

            arc = Arc(Element(center=self.center, radius=self.radius,
                              start_angle=alpha4*180/np.pi,
                              end_angle=alpha1*180/np.pi))
            arc.copy_attributes(self)
            split_arcs.append(arc)
            return split_arcs

        assert(len(points) == 0)
        return []

    def cut_into_halves(self):
        """ return two arcs
        """
        a1 = Arc(Element(center=self.center,
                         radius=self.radius,
                         start_angle=0.0,
                         end_angle=180.0))
        a2 = Arc(Element(center=self.center,
                         radius=self.radius,
                         start_angle=180.0,
                         end_angle=0.0))
        return a1, a2

    def get_angle_of_arc(self):
        return np.pi*2.0

    def is_near(self, n):
        if n[0] > self.center[0] + self.radius + 1e-02:
            return False
        if n[0] < self.center[0] - self.radius - 1e-02:
            return False
        if n[1] > self.center[1] + self.radius + 1e-02:
            return False
        if n[1] < self.center[1] - self.radius - 1e-02:
            return False
        return True

    def __str__(self):
        return "Circle c={}, r={}".format(self.center, self.radius)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self == other

    def __hash__(self):
        """ Override the default hash behavior
            (that returns the id or the object)
        """
        return hash(tuple(sorted(self.__dict__.items())))


#############################
#        Arc (Shape)        #
#############################

class Arc(Circle):
    """a counter clockwise segment of a circle with start and end point"""

    def __init__(self, e, lf=1, rf=np.pi/180,
                 color=None, attr=None, linestyle=None,
                 xoff=0.0, yoff=0.0, rotation=0.0):
        super(self.__class__, self).__init__(e, lf,
                                             xoff=xoff, yoff=yoff,
                                             rotation=rotation)
        self.init_attributes(color, attr, linestyle)

        self.startangle = (e.start_angle + rotation) * rf
        self.endangle = (e.end_angle + rotation) * rf
        if self.endangle < self.startangle:
            if self.endangle < 0:
                self.endangle += 2*np.pi
            elif self.startangle < 0:
                self.startangle += 2*np.pi
            else:
                self.endangle -= 2*np.pi

        self.p1 = (self.center[0] + e.radius*np.cos(self.startangle),
                   self.center[1] + e.radius*np.sin(self.startangle))
        self.p2 = (self.center[0] + e.radius*np.cos(self.endangle),
                   self.center[1] + e.radius*np.sin(self.endangle))
        self.n1 = None
        self.n2 = None

        if hasattr(e, 'rtheta'):
            self.width = e.width
            self.height = e.height
            self.rtheta = e.rtheta
            self.start_param = e.start_param
            self.end_param = e.end_param
        else:
            self.rtheta = None

    def classname(self):
        return "Arc"

    def clone(self):
        alpha_start = alpha_line(self.center, self.p1)
        alpha_end = alpha_line(self.center, self.p2)
        return Arc(Element(center=self.center,
                           radius=self.radius,
                           start_angle=alpha_start*180/np.pi,
                           end_angle=alpha_end*180/np.pi))

    def render(self, renderer, color='blue', with_nodes=False):
        tmp_color = self.get_my_color()
        if not tmp_color:
            tmp_color = color
        tmp_linestyle = self.get_my_linestyle()

        if self.rtheta is None:
            renderer.arc(self.startangle, self.endangle,
                         self.center, self.radius,
                         color=tmp_color,
                         linestyle=tmp_linestyle)
        else:
            renderer.ellipse(self.center, self.width, self.height,
                             self.rtheta,
                             self.start_param,
                             self.end_param,
                             color=tmp_color,
                             linestyle=tmp_linestyle)
        if with_nodes:
            renderer.point(self.p1, 'ro', color)
            renderer.point(self.p2, 'ro', color)

    def center_of_connection(self, ndec=6):
        midangle = middle_angle(self.startangle, self.endangle)
        x, y = self(midangle)
        return (round(x, ndec), round(y, ndec))

    def maxdist(self, r):
        a = np.array([np.linalg.norm(self.center),
                      np.linalg.norm(self.p1),
                      np.linalg.norm(self.p2)])
        return np.max(np.abs(a-r))

    def length(self):
        """returns length of this arc"""
        d = alpha_angle(self.startangle, self.endangle)
        if d > 2*np.pi:
            d -= 2*np.pi
        return self.radius*abs(d)

    def get_alpha(self, n):
        a = alpha_line(n, self.center)
        if points_are_close(n, self.n1):
            alpha1 = normalise_angle(a + np.pi * 0.5)
        elif points_are_close(n, self.n2):
            alpha1 = normalise_angle(a - np.pi * 0.5)
        else:
            alpha1 = 0.0
        return alpha1

    def range(self, step=1.0):
        """returns evenly spaced values"""
        num = max(self.length()/step, 1)
        astep = self.length()/num/self.radius
        s = self.startangle
        d = self.endangle - s
        if d > 2*np.pi:
            d -= 2*np.pi
        alpha = np.arange(s, s+d, astep)
        return self(alpha)

    def __call__(self, alpha):
        """returns x,y coordinates of angle"""
        return (self.center[0] + self.radius*np.cos(alpha),
                self.center[1] + self.radius*np.sin(alpha))

    def overlapping_shape(self, e, rtol=1e-03, atol=1e-03):
        if not isinstance(e, Arc):
            if isinstance(e, Circle):
                return e.overlapping_shape(self, rtol, atol)
            return None

        if not (points_are_close(self.center, e.center) and
                np.isclose(self.radius, e.radius)):
            # Arc (radius/center are not equal)
            return None

        points = []
        if self.is_point_inside(e.p1, rtol=rtol, atol=atol):
            if self.is_point_inside(e.p2, rtol=rtol, atol=atol):
                if e.is_point_inside(self.p2, rtol=rtol, atol=atol):
                    # a Circle
                    points.append(e.p1)
                    points.append(self.p2)
                    points.append(self.p1)
                    points.append(e.p2)
                    points.append(e.p1)
                else:
                    points.append(self.p1)
                    points.append(e.p1)
                    points.append(e.p2)
                    points.append(self.p2)
            else:
                points.append(self.p1)
                points.append(e.p1)
                points.append(self.p2)
                if e.is_point_inside(self.p2, rtol=rtol, atol=atol):
                    points.append(e.p2)

        elif self.is_point_inside(e.p2, rtol=rtol, atol=atol):
            points.append(e.p1)
            points.append(self.p1)
            points.append(e.p2)
            points.append(self.p2)

        elif e.is_point_inside(self.p1, rtol=rtol, atol=atol):
            points.append(e.p1)
            points.append(self.p1)
            points.append(self.p2)
            if e.is_point_inside(self.p2, rtol=rtol, atol=atol):
                points.append(e.p2)

        elif e.is_point_inside(self.p2, rtol=rtol, atol=atol):
            if not points_are_close(self.p1, e.p1, rtol=rtol, atol=atol):
                logger.error("FATAL ERROR in overlapping_shape() of Arc")

                raise ValueError('FATAL ERROR in overlapping_shape() of Arc')
            points.append(e.p1)
            points.append(self.p2)
            points.append(e.p2)

        elif points_are_close(self.p1, e.p1) and points_are_close(self.p2, e.p2):
            # self and e are identical
            return [e]

        else:
            # no overlapping elements
            return None

        return self.create_arcs(points)

    def intersect_line(self, line, rtol=1e-03, atol=1e-03, include_end=False):
        """ Von einem Arc-Objekt und einem Line-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben
        """
        points = super(Arc, self).intersect_line(line, rtol, atol, include_end)
        if not points:
            return []

        # all possible points have been found
        # Lets see if they are on a arc
        remaining_points = []
        for p in points:
            if self.is_point_inside(p,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=include_end):
                remaining_points.append(p)
        return remaining_points

    def intersect_arc(self, arc, rtol=1e-03, atol=1e-03, include_end=False):
        """ Von zwei Arc-Objekten werden die Schnittpunkte bestimmt und in
            einer Liste ausgegeben.
        """
        assert(isinstance(arc, Arc))
        points = self.intersect_circle(arc, rtol, atol, include_end)

        # Check if the points are on a arc
        # (has been assumed as a circle)
        remaining_points = []
        for p in points:
            if arc.is_point_inside(p,
                                   rtol=rtol,
                                   atol=atol,
                                   include_end=include_end):
                remaining_points.append(p)
        return remaining_points

    def intersect_circle(self, circle, rtol=1e-03,
                         atol=1e-03, include_end=False):
        """ return the list of intersection points
        """
        if points_are_close(self.center, circle.center, rtol, atol):
            if np.isclose(self.radius, circle.radius):
                if include_end:
                    return [self.p1, self.p2]
            # no intersection with different radius but equal center
            return []

        points = super(Arc, self).intersect_circle(
            circle, rtol, atol, include_end)

        # Intersection points exist. Take the ones on the arc
        remaining_points = []
        for p in points:
            if self.is_point_inside(p,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=include_end):
                remaining_points.append(p)
        return remaining_points

    def split(self, points, rtol=1e-03, atol=1e-03, mdec=0):
        """ return a list of arcs by splitting
        """
        points_inside = [p
                         for p in points
                         if self.is_point_inside(p,
                                                 rtol=rtol,
                                                 atol=atol)]
        if len(points_inside) == 1:
            p = points_inside[0]
            split_arcs = []
            alpha = alpha_line(self.center, p)
            arc = Arc(Element(center=self.center, radius=self.radius,
                              start_angle=self.startangle*180/np.pi,
                              end_angle=alpha*180/np.pi))
            arc.copy_attributes(self)
            split_arcs.append(arc)

            arc = Arc(Element(center=self.center, radius=self.radius,
                              start_angle=alpha*180/np.pi,
                              end_angle=self.endangle*180/np.pi))
            arc.copy_attributes(self)
            split_arcs.append(arc)
            return split_arcs

        assert(len(points_inside) == 0)
        return []

    def cut_into_halves(self):
        """ return two arcs
        """
        midangle = middle_angle(self.startangle, self.endangle)
        a1 = Arc(Element(center=self.center,
                         radius=self.radius,
                         start_angle=self.startangle*180/np.pi,
                         end_angle=midangle*180/np.pi))
        a2 = Arc(Element(center=self.center,
                         radius=self.radius,
                         start_angle=midangle*180/np.pi,
                         end_angle=self.endangle*180/np.pi))
        return a1, a2

    def concatenate_arc(self, n1, n2, el,
                        rtol=1e-03, atol=1e-03,
                        overlapping=False):
        if not points_are_close(self.center, el.center):
            return None
        if not np.isclose(self.radius, el.radius):
            return None

        my_start = normalise_angle(self.startangle)
        my_end = normalise_angle(self.endangle)
        el_start = normalise_angle(el.startangle)
        el_end = normalise_angle(el.endangle)

        logger.debug("begin of concatenate_arc")
        logger.debug("my startangle: %s,  endangle: %s", my_start, my_end)
        logger.debug("el startangle: %s,  endangle: %s", el_start, el_end)

        start_angle = None
        end_angle = None

        if np.isclose(my_start, el_end, rtol=rtol, atol=atol):
            if is_angle_inside(my_start, my_end, el_start):  # Circle
                logger.debug("1: new startangle is 0 (Circle)")
                start_angle = 0.0
                end_angle = 0.0
            else:
                logger.debug("1: new startangle is el: %s", el_start)
                start_angle = el_start
                end_angle = my_end
        elif np.isclose(el_start, my_end, rtol=rtol, atol=atol):
            logger.debug("2: new startangle is my: %s", my_start)
            start_angle = my_start
            end_angle = el_end

        if start_angle is None:
            if overlapping:
                if is_angle_inside(my_start, my_end, el_start):
                    logger.debug("3: new startangle is my: %s", my_start)
                    start_angle = my_start
                elif is_angle_inside(el_start, el_end, my_start):
                    logger.debug("4: new startangle is el: %s", el_start)
                    start_angle = el_start
                if start_angle is None:
                    return None

                if is_angle_inside(my_start, my_end, el_end):
                    end_angle = my_end
                elif is_angle_inside(el_start, el_end, my_end):
                    end_angle = el_end
                if end_angle is None:
                    return None
            else:
                return None

        logger.debug("new startangle: %s,  endangle: %s", start_angle, end_angle)

        if np.isclose(start_angle, end_angle):
            # it's a circle
            return Circle(Element(center=self.center, radius=self.radius))

        logger.debug("end of concatenate_arc")
        return Arc(Element(center=self.center,
                           radius=self.radius,
                           start_angle=start_angle,
                           end_angle=end_angle),
                   rf=1)

    def concatenate_circle(self, n1, n2, el,
                           rtol=1e-03, atol=1e-03,
                           overlapping=False):
        if not points_are_close(self.center, el.center):
            return None
        if not np.isclose(self.radius, el.radius):
            return None
        # it's a circle
        return Circle(Element(center=self.center, radius=self.radius))

    def is_point_inside(self, p,
                        rtol=1e-03,
                        atol=1e-03,
                        include_end=False,
                        ignore_end=False,
                        mdec=0):
        """ returns true if p is on arc
        """
        # logger.debug("is_point_inside: p=%s", p)
        d = distance(p, self.center)
        if not np.isclose(d, self.radius, rtol=rtol, atol=atol):
            # logger.debug(" <== RADIUS %s, DISTANCE %s",
            #              self.radius, d)
            return False
        if points_are_close(p, self.p1, rtol=rtol, atol=atol):
            # logger.debug(" <== CLOSE TO P1 %s: rtol=%s, atol=%s",
            #              self.p1, rtol, atol)
            return include_end
        elif points_are_close(p, self.p2, rtol=rtol, atol=atol):
            # logger.debug(" <== CLOSE TO P2 %s: rtol=%s, atol=%s",
            #              self.p2, rtol, atol)
            return include_end
        elif points_are_close(self.p1, self.p2, rtol=rtol, atol=atol):
            # logger.debug(" <== P1 AND P2 CLOSE TOGETHER")
            return False

        alpha_p1 = alpha_line(self.center, self.p1)
        alpha_p2 = alpha_line(self.center, self.p2)
        alpha_p = alpha_line(self.center, p)
        alpha_inside = is_angle_inside(alpha_p1, alpha_p2, alpha_p)
        # logger.debug("is_point_inside: %s (%s, %s ,%s)",
        #              alpha_inside, alpha_p1, alpha_p2, alpha_p)
        return alpha_inside

    def is_angle_inside(self, alpha, rtol=1e-03, atol=1e-03,
                        include_end=False):
        """ returns True if alpha is between start and end angle
        """
        return is_angle_inside(self.startangle, self.endangle, alpha)

    def transform(self, T, alpha, ndec):
        super(Arc, self).transform(T, alpha, ndec)
        p1, p2 = ((self.p1[0]-self.center[0],
                   self.p1[1]-self.center[1]),
                  (self.p2[0]-self.center[0],
                   self.p2[1]-self.center[1]))

        self.startangle = np.arctan2(p1[1], p1[0])
        self.endangle = np.arctan2(p2[1], p2[0])
        if self.rtheta is not None:
            self.rtheta = self.rtheta + alpha
        return self

    def correct(self, src_alpha, dest_alpha, ndec):
        super(Arc, self).correct(src_alpha, dest_alpha, ndec)

        p1, p2 = ((self.p1[0]-self.center[0],
                   self.p1[1]-self.center[1]),
                  (self.p2[0]-self.center[0],
                   self.p2[1]-self.center[1]))

        self.startangle = np.arctan2(p1[1], p1[0])
        self.endangle = np.arctan2(p2[1], p2[0])
        if self.rtheta is not None:
            self.rtheta = self.rtheta + dest_alpha
        return self

    def minmax(self):
        """ Die Funktion bestimmt das Minimum und Maximum auf der x- und der
            y-Achse (return [<min-x>, <max-x>, <min-y>, <max-y>])
        """
        mm = [min(self.p1[0], self.p2[0]), max(self.p1[0], self.p2[0]),
              min(self.p1[1], self.p2[1]), max(self.p1[1], self.p2[1])]

        p = [self.center[0]-self.radius, self.center[1]]
        if p[0] < mm[0]:
            a = alpha_line(self.center, p)
            if self.is_angle_inside(a, 0.00001):
                mm[0] = p[0]

        p = [self.center[0]+self.radius, self.center[1]]
        if p[0] > mm[1]:
            a = alpha_line(self.center, p)
            if self.is_angle_inside(a, 0.00001):
                mm[1] = p[0]

        p = [self.center[0], self.center[1]-self.radius]
        if p[1] < mm[2]:
            a = alpha_line(self.center, p)
            if self.is_angle_inside(a, 0.00001):
                mm[2] = p[1]

        p = [self.center[0], self.center[1]+self.radius]
        if p[1] > mm[3]:
            a = alpha_line(self.center, p)
            if self.is_angle_inside(a, 0.00001):
                mm[3] = p[1]
        return mm

    def minmax_from_center(self, center):
        """ Die Funktion ermittelt den minimalen und maximalen
            Abstand vom Center
        """
        d = distance(center, self.center)
        if np.isclose(d, 0.0):
            return (self.radius, self.radius)

        angle = alpha_line(center, self.center)
        dist_min = abs(d - self.radius)
        dist_max = d + self.radius

        pmax = point(center, d + self.radius, angle)
        alpha_pmax = alpha_line(self.center, pmax)
        if not self.is_angle_inside(alpha_pmax, 1e-08):
            dist_max = max(distance(center, self.p1),
                           distance(center, self.p2))

        pmin = point(center, d - self.radius, angle)
        alpha_pmin = alpha_line(self.center, pmin)

        if not self.is_angle_inside(alpha_pmin, 1e-08):
            dist_min = min(distance(center, self.p1),
                           distance(center, self.p2))

        return (dist_min, dist_max)

    def minmax_angle_from_center(self, center):
        d = distance(center, self.center)
        r = self.radius
        v = d**2 - r**2
        if less_equal(v, 0.0):
            points = []
        else:
            r2 = np.sqrt(v)
            circ = Circle(Element(center=center, radius=r2))
            points = self.intersect_circle(circ, include_end=True)

        points.append(self.p2)

        alpha_min = alpha_line(center, self.p1)
        alpha_max = alpha_min

        for p in points:
            alpha_p = alpha_line(center, p)
            alpha_min = min_angle(alpha_min, alpha_p)
            alpha_max = max_angle(alpha_max, alpha_p)

        return (alpha_min, alpha_max)

    def get_nodes(self, parts=8, render=False):
        """ Die Funktion liefert eine Liste von virtuellen Nodes, welche man
            zum Rechnen der convex_hull() benötigt.
        """
        if render and self.rtheta is not None:
            theta = np.arange(self.start_param,
                              self.end_param,
                              0.1)
            x = 0.5 * self.width * np.cos(theta)
            y = 0.5 * self.height * np.sin(theta)
            R = np.array([
                [np.cos(self.rtheta), -np.sin(self.rtheta)],
                [np.sin(self.rtheta),  np.cos(self.rtheta)]])
            x, y = np.dot(R, np.array([x, y]))
            x += self.center[0]
            y += self.center[1]
            nodes = list(zip(x, y))
            nodes.append(self.p2)
            return nodes

        return (p for p in points_on_arc(self.center, self.radius,
                                         self.startangle,
                                         self.endangle,
                                         parts=parts))

    def get_angle_of_arc(self):
        return get_angle_of_arc(self.startangle, self.endangle)

    def __str__(self):
        return "Arc c={},\n\t r={},\n\t start={},\n\t end={},\n\t p1={},\n\t p2={},\n\t n1={},\n\t n2={}".\
            format(self.center,
                   self.radius, self.startangle,
                   self.endangle,
                   self.p1,
                   self.p2,
                   self.n1,
                   self.n2)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self == other

    def __hash__(self):
        """ Override the default hash behavior
            (that returns the id or the object)
        """
        return hash(tuple(sorted(self.__dict__.items())))


#############################
#        Line (Shape)       #
#############################

class Line(Shape):
    """straight connection between start and end point"""

    def __init__(self, e, lf=1,
                 color=None, attr=None, linestyle=None,
                 xoff=0.0, yoff=0.0, rotation=0.0):
        self.init_attributes(color, attr, linestyle)
        if rotation != 0.0:
            alpha = rotation * np.pi/180
            T = np.array(((np.cos(alpha), -np.sin(alpha)),
                          (np.sin(alpha), np.cos(alpha))))
            start = self.rotate(T, e.start)
            end = self.rotate(T, e.end)
        else:
            start = e.start
            end = e.end
        self.p1 = lf*start[0] + xoff, lf*start[1] + yoff
        self.p2 = lf*end[0] + xoff, lf*end[1] + yoff
        self.n1 = None
        self.n2 = None

    def classname(self):
        return "Line"

    def clone(self):
        return Line(Element(start=self.p1,
                            end=self.p2))

    def render(self, renderer, color='blue', with_nodes=False):
        tmp_color = self.get_my_color()
        if not tmp_color:
            tmp_color = color
        tmp_linestyle = self.get_my_linestyle()

        renderer.line(self.p1, self.p2,
                      color=tmp_color,
                      linestyle=tmp_linestyle,
                      e=self)
        if with_nodes:
            renderer.point(self.p1, 'ro', tmp_color)
            renderer.point(self.p2, 'ro', tmp_color)

    def center_of_connection(self, ndec=6):
        x = (self.p1[0]+self.p2[0])/2
        y = (self.p1[1]+self.p2[1])/2
        return (x, y)

    def maxdist(self, r):
        r1 = np.linalg.norm(self.p1)
        r2 = np.linalg.norm(self.p2)
        return max([abs(r1-r), abs(r2-r)])

    def length(self):
        return np.sqrt(self.dx()**2 + self.dy()**2)

    def get_alpha(self, n):
        px = self.center_of_connection()
        return alpha_line(px, n)

    def get_positive_angle(self):
        return elevation_angle(alpha_line(self.p1, self.p2))

    def range(self, step=1.0):
        """returns evenly spaced values"""
        num = max(int(self.length()/step), 1)
        if np.isclose(self.dx(), 0):
            if np.isclose(self.dy(), 0):
                return ((self.xmin(),), (self.ymin(),))
            y = np.arange(self.ymin(), self.ymax(), self.length()/num)
            return np.array((self.xmin()*np.ones(len(y)), y))
        xstep = np.sqrt((self.length()/num)**2/(1 + self.dy()**2/self.dx()**2))
        x = np.arange(self.xmin(), self.xmax(), xstep)
        return (x, self(x))

    def __call__(self, x):
        """returns y coordinate of x"""
        if np.isclose(self.dx(), 0):
            return float('nan')
        return self.dy()/self.dx()*(x - self.p1[0]) + self.p1[1]

    def overlapping_shapes(self, n, e, rtol=1e-03, atol=1e-03):
        if not isinstance(e, Line):
            return False

        if nodes_are_equal(n, self.n1):
            my_m = alpha_line(n, self.n2)
        else:
            my_m = alpha_line(n, self.n1)

        if nodes_are_equal(n, e.n1):
            other_m = alpha_line(n, e.n2)
        else:
            other_m = alpha_line(n, e.n1)

        return np.isclose(my_m, other_m, rtol=rtol, atol=atol)

    def intersect_line(self, line, rtol=1e-03, atol=1e-03,
                       include_end=False, all=False):
        """ Von zwei Line-Objekten wird der Schnittpunkt bestimmt und in
            einer Liste ausgegeben.
        """
        point = []
        m_L1 = self.m()
        m_L2 = line.m()
        if m_L1 is None:
            if m_L2 is None:
                return []
            else:
                y = line_n([line.p1[0]-self.p1[0], line.p1[1]], m_L2)
                point = (self.p1[0], y)
        else:
            if m_L2 is None:
                y = line_n([self.p1[0]-line.p1[0], self.p1[1]], m_L1)
                point = (line.p1[0], y)
            else:
                if np.isclose(m_L1, m_L2):
                    return []
                else:
                    point = lines_intersect_point(self.p1, m_L1, self.n(m_L1),
                                                  line.p1, m_L2, line.n(m_L2))

        if all:
            return[point]

        if line.is_point_inside(point,
                                rtol=rtol,
                                atol=atol,
                                include_end=include_end):
            if self.is_point_inside(point,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=include_end):
                return [point]
        return []

    def intersect_arc(self, arc, rtol=1e-03, atol=1e-03, include_end=False):
        """ Von einem Line-Objekt und einem Arc-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben.
        """
        return arc.intersect_line(self, rtol, atol, include_end)

    def intersect_circle(self, circle, rtol=1e-03, atol=1e-03,
                         include_end=False):
        """ Von einem Line-Objekt und einem Circle-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben.
        """
        return circle.intersect_line(self, rtol, atol, include_end)

    def split(self, points, rtol=1e-03, atol=1e-03, mdec=0):
        """ Die Funktion splittet das Line-Objekt an den vorgegebenen Punkten
            und gibt eine Liste der neu enstandenen Elemente aus.
        """
        points_inside = [(distance(p, self.p1), p)
                         for p in points if self.is_point_inside(p,
                                                                 rtol=rtol,
                                                                 atol=atol,
                                                                 mdec=mdec)]
        if len(points_inside) > 0:
            points_inside.append((0.0, self.p1))
            points_inside.append((distance(self.p1, self.p2), self.p2))
            points_inside.sort()

            split_lines = []
            p_start = None
            for d, p in points_inside:
                if p_start is not None:
                    line = Line(Element(start=p_start, end=p))
                    line.copy_attributes(self)
                    split_lines.append(line)

                p_start = p
            return split_lines
        return []

    def cut_into_halves(self):
        """ return two lines
        """
        pm = middle_point_of_line(self.p1, self.p2)
        l1 = Line(Element(start=self.p1, end=pm))
        l2 = Line(Element(start=pm, end=self.p2))
        return l1, l2

    def concatenate_line(self, n1, n2, el,
                         rtol=1e-05, atol=1e-05,
                         mdec=0,
                         overlapping=False):
        none_val = 9999999.0
        my_m = self.m(none_val, dec=mdec)
        el_m = el.m(none_val, dec=mdec)
        #logger.debug("concatenate_line: my m=%s, el m=%s", my_m, el_m)
        if not np.isclose(my_m, el_m):
            #logger.debug("concatenate_line #1: m %s and %s are not close", my_m, el_m)
            return None
        if n1 and n2:
            #logger.debug("concatenate_line #2: Nodes %s and %s", n1, n2)
            return Line(Element(start=n1, end=n2))

        my_p1, my_p2 = self.points_sorted(rtol=rtol, atol=atol)
        el_p1, el_p2 = el.points_sorted(rtol=rtol, atol=atol)

        logger.debug("my p1: %s,  p2: %s", my_p1,my_p2)
        logger.debug("el p1: %s,  p2: %s", el_p1,el_p2)

        if points_are_close(my_p2, el_p1, rtol=rtol, atol=atol):
            return Line(Element(start=my_p1, end=el_p2))
        if points_are_close(my_p1, el_p2, rtol=rtol, atol=atol):
            return Line(Element(start=el_p1, end=my_p2))

        if overlapping:
            if points_are_close(my_p1, el_p1, rtol=rtol, atol=atol):
                my_d = distance(my_p1, my_p2)
                el_d = distance(el_p1, el_p2)
                if el_d > my_d:
                    return Line(Element(start=el_p1, end=el_p2))
                else:
                    return Line(Element(start=my_p1, end=my_p2))

            if points_are_close(my_p2, el_p2, rtol=rtol, atol=atol):
                my_d = distance(my_p1, my_p2)
                el_d = distance(el_p1, el_p2)
                if el_d > my_d:
                    return Line(Element(start=el_p1, end=el_p2))
                else:
                    return Line(Element(start=my_p1, end=my_p2))

            start = None
            if self.is_point_inside(el_p1,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=True,
                                    mdec=mdec):
                start = my_p1
            elif el.is_point_inside(my_p1,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=True,
                                    mdec=mdec):
                start = el_p1
            if start is None:
                logger.debug("concatenate_line #4: no start")
                return None
            end = None
            if self.is_point_inside(el_p2,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=True,
                                    mdec=mdec):
                end = my_p2
            elif el.is_point_inside(my_p2,
                                    rtol=rtol,
                                    atol=atol,
                                    include_end=True,
                                    mdec=mdec):
                end = el_p2
            if end is None:
                logger.debug("concatenate_line #4: no end")
                return None

            return Line(Element(start=start, end=end))

        return None

    def is_point_inside(self, point,
                        rtol=1e-03,
                        atol=1e-03,
                        include_end=False,
                        ignore_end=False,
                        mdec=0):
        """ returns True if point is between start and end point of the line
        """
        if not ignore_end:
            if points_are_close(point, self.p1, rtol, atol):
                return include_end
            if points_are_close(point, self.p2, rtol, atol):
                return include_end

        elevation_tol = np.pi / 360 / 8  # 1/8 degree
        elevation_line = elevation_angle(alpha_line(self.p1, self.p2))
        elevation_point1 = elevation_angle(alpha_line(self.p1, point))
        elevation_point2 = elevation_angle(alpha_line(point, self.p2))

        if not (np.isclose(elevation_line, elevation_point1, rtol=rtol, atol=atol) or \
                np.isclose(elevation_line, elevation_point2, rtol=rtol, atol=atol)):
            return False

        length = distance(self.p1, self.p2)
        dist_p1 = distance(self.p1, point)
        dist_p2 = distance(self.p2, point)
        if dist_p1 > length or dist_p2 > length:
            return False
        return True

    def minmax(self):
        """ Die Funktion bestimmt das Minimum und Maximum auf der x- und der
            y-Achse (return [<min-x>, <max-x>, <min-y>, <max-y>])
        """
        return [min(self.p1[0], self.p2[0]), max(self.p1[0], self.p2[0]),
                min(self.p1[1], self.p2[1]), max(self.p1[1], self.p2[1])]

    def minmax_from_center(self, center):
        """ Die Funktion ermittelt den minimalen und maximalen Abstand vom Center
        """
        dist_max = max(distance(center, self.p1), distance(center, self.p2))
        dist_min = min(distance(center, self.p1), distance(center, self.p2))
        m = line_m(self.p1, self.p2)
        n = line_n(self.p1, m)
        p = intersect_point(center, self.p1, m, n)

        if self.is_point_inside(p,
                                rtol=1e-03,
                                atol=1e-03):
            dist_min = min(distance(center, p), dist_min)

        return (dist_min, dist_max)

    def minmax_angle_from_center(self, center):
        if points_are_close(center, self.p1):
            alpha_p1 = None
        else:
            alpha_p1 = alpha_line(center, self.p1)
        if points_are_close(center, self.p2):
            alpha_p2 = None
        else:
            alpha_p2 = alpha_line(center, self.p2)
        if alpha_p1 is None:
            assert(alpha_p2 is not None)
            return (alpha_p2, alpha_p2)
        if alpha_p2 is None:
            return (alpha_p1, alpha_p1)
        if alpha_angle(alpha_p1, alpha_p2) < np.pi:
            return (alpha_p1, alpha_p2)
        else:
            return (alpha_p2, alpha_p1)

    def get_nodes(self, parts=8, render=False):
        """ Die Funktion liefert eine Liste von virtuellen Nodes, welche man
            zum Rechnen der convex_hull() benötigt.
        """
        return (self.p1, self.p2)

    def get_angle_of_arc(self):
        return 0.0

    def is_near(self, n):
        if n[0] > max(self.p1[0], self.p2[0]) + 1e-02:
            return False
        if n[0] < min(self.p1[0], self.p2[0]) - 1e-02:
            return False
        if n[1] > max(self.p1[1], self.p2[1]) + 1e-02:
            return False
        if n[1] < min(self.p1[1], self.p2[1]) - 1e-02:
            return False
        return True

    def __str__(self):
        return "Line p1={}, p2={}, n1={}, n2={}".format(self.p1,
                                                        self.p2,
                                                        self.n1,
                                                        self.n2)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self == other

    def __hash__(self):
        """ Override the default hash behavior
            (that returns the id or the object)
        """
        return hash(tuple(sorted(self.__dict__.items())))


#############################
#       Point (Shape)       #
#############################

class Point(Shape):
    """ used for plotting only """

    def __init__(self, p):
        self.p1 = p

    def classname(self):
        return "Point"

    def render(self, renderer):
        renderer.point(self.p1)

    def transform(self, T, alpha, ndec):
        n = T.dot(np.array((self.p1[0], self.p1[1])))
        self.p1 = (n[0], n[1])
        return self


def is_Circle(e):
    return isinstance(e, Circle) and not isinstance(e, Arc)


def is_Arc(e):
    return isinstance(e, Arc)


def is_Line(e):
    return isinstance(e, Line)
