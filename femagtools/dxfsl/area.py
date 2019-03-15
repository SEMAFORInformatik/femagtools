# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.area
  ~~~~~~~~~~~~~~~~~~~~~

  areas are regions surrounded by a collection of shapes

  Authors: Ronald Tanner, beat Holm
"""
from __future__ import print_function
import numpy as np
import networkx as nx
import logging
from .functions import less_equal, less, greater_equal, greater
from .functions import distance, alpha_angle, alpha_line, min_angle, max_angle
from .functions import point, line_m, line_n, intersect_point, points_are_close
from .functions import middle_angle, part_of_circle, is_same_angle
from .shape import Element, Shape, Line, Arc, Circle

logger = logging.getLogger('femagtools.area')


#############################
#            Area           #
#############################

area_number = 0


class Area(object):
    def __init__(self, area, center, sym_tolerance):
        self.area = area
        self.type = 0  # material
        self.phi = 0.0
        self.min_angle = 0.0
        self.max_angle = 0.0
        self.min_air_angle = 0.0
        self.max_air_angle = 0.0
        self.close_to_ag = False
        self.close_to_startangle = False
        self.close_to_endangle = False
        self.mag_rectangle = False
        self.min_dist = 99999.0
        self.max_dist = 0.0
        self.height = 0.0
        self.alpha = 0.0
        self.count = 1
        self.equal_areas = []
        self.delta = 0.0
        self.start = 0.0
        self.sym_startangle = 0.0
        self.sym_endangle = 0.0
        self.sym_type = 0
        self.symmetry = 0
        self.sym_tolerance = sym_tolerance
        self.calc_signature(center)
        global area_number
        area_number += 1
        self.id = area_number

    def identifier(self):
        return "{}-{}".format(self.id, self.type)

    def number_of_elements(self):
        return len(self.area)

    def elements(self):
        return self.area

    def nodes(self):
        if len(self.area) == 0:
            return

        for e in self.area:
            yield e.p1
            yield e.p2

    def virtual_nodes(self):
        if len(self.area) < 2:
            return

        prev_nodes = [n for n in self.area[0].get_nodes(parts=64)]
        next_nodes = [n for n in self.area[1].get_nodes(parts=64)]
        if points_are_close(prev_nodes[0], next_nodes[0], 1e-03, 1e-01):
            prev_nodes = prev_nodes[::-1]
        elif points_are_close(prev_nodes[0], next_nodes[-1], 1e-03, 1e-01):
            prev_nodes = prev_nodes[::-1]
            next_nodes = next_nodes[::-1]
        elif points_are_close(prev_nodes[-1], next_nodes[-1], 1e-03, 1e-01):
            next_nodes = next_nodes[::-1]
        elif not points_are_close(prev_nodes[-1], next_nodes[0], 1e-03, 1e-01):
            assert(False)
        last_point = next_nodes[-1]
        for n in prev_nodes:
            yield n
        for n in next_nodes:
            yield n

        for e in self.area[2::]:
            next_nodes = [n for n in e.get_nodes(parts=64)]

            if points_are_close(next_nodes[-1], last_point, 1e-03, 1e-01):
                next_nodes = next_nodes[::-1]
            for n in next_nodes:
                yield n
            last_point = next_nodes[-1]

    def name(self):
        if self.type == 1:
            return 'Iron'
        if self.type == 2:
            return 'windings'
        if self.type == 3 or self.type == 4:
            return 'magnet'
        if self.type == 5:
            return 'StJo'
        if self.type == 6:
            return 'StZa'
        return ''

    def color(self):
        if self.type == 1:
            return 'blue'
        if self.type == 2:
            return 'green'
        if self.type == 3 or self.type == 4:
            return 'red'
        if self.type == 5:
            return 'cyan'
        if self.type == 6:
            return 'blue'
        return 'white'

    def is_iron(self):
        return self.type == 1 or self.type == 5 or self.type == 6

    def is_stator_iron_yoke(self):
        return self.type == 5

    def is_stator_iron_shaft(self):
        return self.type == 6

    def is_rotor_iron(self):
        return self.type == 1

    def is_winding(self):
        return self.type == 2

    def is_magnet(self):
        return self.type == 3 or self.type == 4

    def is_air(self):
        return self.type == 0

    def set_type(self, t):
        self.type = t

    def calc_signature(self, center):
        if not self.area:
            return

        s = self.area[0]
        mm_angle = s.minmax_angle_from_center(center)
        self.min_angle = mm_angle[0]
        self.max_angle = mm_angle[1]

        for s in self.area:
            mm_dist = s.minmax_from_center(center)
            self.min_dist = min(self.min_dist, mm_dist[0])
            self.max_dist = max(self.max_dist, mm_dist[1])
            self.height = self.max_dist - self.min_dist

            mm_angle = s.minmax_angle_from_center(center)
            self.min_angle = min_angle(self.min_angle, mm_angle[0])
            self.max_angle = max_angle(self.max_angle, mm_angle[1])

        self.alpha = round(alpha_angle(self.min_angle, self.max_angle), 3)

    def minmax_angle_dist_from_center(self, center, dist):
        circ = Circle(Element(center=center, radius=dist))
        s = self.area[0]
        my_min_angle = self.max_angle
        my_max_angle = self.min_angle
        mm_angle = None
        for s in self.area:
            mm_angle = s.minmax_angle_dist_from_center(my_min_angle,
                                                       my_max_angle,
                                                       center, circ)
            if mm_angle:
                my_min_angle = min_angle(my_min_angle, mm_angle[0])
                my_max_angle = max_angle(my_max_angle, mm_angle[1])
        return (my_min_angle, my_max_angle)

    def is_inside(self, area):
        if less_equal(area.min_dist, self.min_dist):
            return False
        if greater_equal(area.max_dist, self.max_dist):
            return False
        if less_equal(area.min_angle, self.min_angle):
            return False
        if greater_equal(area.max_angle, self.max_angle):
            return False
        return True

    def is_touching(self, area):
        for n in self.nodes():
            x = [p for p in area.nodes() if points_are_close(n, p)]
            if x:
                return True
        return False

    def has_connection(self, geom, a, ndec):
        assert(self.area)
        assert(a.area)
        n1 = self.area[0].node1(ndec)
        if not geom.g.has_node(n1):
            n = geom.find_nodes(n1)
            if not n:
                logger.warn("FATAL: node {} not available".format(n1))
                return False
            n1 = n[0]

        n2 = a.area[0].node2(ndec)
        if not geom.g.has_node(n2):
            n = geom.find_nodes(n2)
            if not n:
                logger.warn("FATAL: node {} not available".format(n2))
                return False
            n2 = n[0]

        try:
            return nx.has_path(geom.g, n1, n2)
        except nx.NetworkXError:
            logger.warn("has_path() failed")
            return False

    def get_lowest_gap_list(self, a, center, radius, rightangle, leftangle):
        gap_list = []
        for p1 in self.nodes():
            for p2 in a.nodes():
                d = distance(p1, p2)
                gap_list.append((d, (p1, p2)))

        d, p1, p2 = a.get_nearest_point(center, radius, rightangle)
        gap_list.append((d, (p1, p2)))
        d, p1, p2 = a.get_nearest_point(center, radius, leftangle)
        gap_list.append((d, (p1, p2)))
        return gap_list

    def get_nearest_point(self, center, radius, angle):
        axis_p = point(center, radius, angle)
        axis_m = line_m(center, axis_p)
        axis_n = line_n(center, axis_m)

        the_area_p = None
        the_axis_p = None
        dist = 99999
        for n in self.nodes():
            p = intersect_point(n, center, axis_m, axis_n)
            d = distance(n, p)
            if d < dist:
                dist = d
                the_area_p = n
                the_axis_p = p

        return (dist,
                (the_axis_p[0], the_axis_p[1]),
                (the_area_p[0], the_area_p[1]))

    def is_equal(self, a, sym_tolerance):
        if sym_tolerance > 0.0:
            if np.isclose(round(self.min_dist, 4),
                          round(a.min_dist, 4),
                          1e-03, sym_tolerance) and \
               np.isclose(round(self.max_dist, 4),
                          round(a.max_dist, 4),
                          1e-03, sym_tolerance) and \
               np.isclose(round(self.alpha, 3),
                          round(a.alpha, 3),
                          1e-02, 0.01):
                return True
        else:
            if np.isclose(round(self.min_dist, 2),
                          round(a.min_dist, 2)) and \
               np.isclose(round(self.max_dist, 2),
                          round(a.max_dist, 2)) and \
               np.isclose(round(self.alpha, 3),
                          round(a.alpha, 3), 1e-02, 0.001):
                return True
        return False

    def is_identical(self, area):
        if np.isclose(self.min_dist, area.min_dist) and \
           np.isclose(self.max_dist, area.max_dist) and \
           np.isclose(self.alpha, area.alpha) and \
           np.isclose(self.min_angle, area.min_angle) and \
           np.isclose(self.max_angle, area.max_angle):
            return True
        return False

    def increment(self, a):
        if self.is_identical(a):
            return

        for area in self.equal_areas:
            if area.is_identical(a):
                return

        self.count += 1
        self.equal_areas.append(a)

    def set_delta(self):
        logger.debug("begin set_delta of {}".format(self.id))
        self.delta = 0.0
        self.symmetry = 0

        if len(self.equal_areas) < 2:
            # Mit zwei Objekten lässt sich das Teil nur noch halbieren. Das
            # wird zum Schluss sowieso versucht.
            logger.debug("end set_delta: zuwenig Gleiche")
            return

        sorted_areas = []
        sorted_areas.append((self.min_angle, self))
        for a in self.equal_areas:
            sorted_areas.append((a.min_angle, a))
        sorted_areas.sort()

        delta = {}
        prev_angle = sorted_areas[0][0]
        for angle, area in sorted_areas[1:]:
            d = round(alpha_angle(prev_angle, angle), 2)
            if d in delta:
                delta[d] += 1
            else:
                delta[d] = 1
            prev_angle = angle

        delta_sorted = list([v, k] for (k, v) in delta.items())
        logger.debug(" - delta: {}".format(delta_sorted))

        if len(delta_sorted) == 1:
            # simple case: all have the same angle
            self.delta = alpha_angle(sorted_areas[0][1].min_angle,
                                     sorted_areas[1][1].min_angle)
            self.start = middle_angle(sorted_areas[0][1].max_angle,
                                      sorted_areas[1][1].min_angle)
            self.sym_type = 3
            self.symmetry = part_of_circle(0.0, self.delta, 1)
            logger.debug("end set_delta: simple case")
            return

        logger.debug("end set_delta: {} deltas, {} areas"
                     .format(len(delta_sorted), len(self.equal_areas)))

        if len(delta_sorted) > 2:
            # Mehr als 2 Winkel untersuchen wir (noch) nicht. Wir brechen
            # die Suche nach dem richtigen Winkel ab.
            logger.debug("end set_delta: zuviele Winkel")
            return

        # Bei 2 verschiedenen Winkeln werden die näher beieinander liegenden
        # Objekte zusammen genommen.

        if len(self.equal_areas) < 4:
            # Wenn nicht mehr als 4 Objekte vorhanden sind, brechen wir auch
            # ab.
            logger.debug("end set_delta: zuwenig areas")
            return

        if np.isclose(delta_sorted[1][1],
                      delta_sorted[0][1]*2, atol=0.01):
            # Lets hope we have unreqognised areas inbetween
            percent = 1.0
        else:
            percent = delta_sorted[0][0] / (len(self.equal_areas)+1)

        if percent > 0.75:
            # lets assume we only have one angle
            self.delta = alpha_angle(sorted_areas[0][1].min_angle,
                                     sorted_areas[1][1].min_angle)
            self.start = middle_angle(sorted_areas[0][1].max_angle,
                                      sorted_areas[1][1].min_angle)
            self.sym_type = 2
            self.symmetry = part_of_circle(0.0, self.delta, 1)
            logger.debug("end set_delta: {} Prozent gleiche deltas"
                         .format(percent))
            return

        # Lets hope the distances are changing
        self.delta = alpha_angle(sorted_areas[0][1].min_angle,
                                 sorted_areas[2][1].min_angle)
        self.sym_type = 1
        self.symmetry = part_of_circle(0.0, self.delta, 1)
        delta_1 = alpha_angle(sorted_areas[0][1].min_angle,
                              sorted_areas[1][1].min_angle)
        delta_2 = alpha_angle(sorted_areas[1][1].min_angle,
                              sorted_areas[2][1].min_angle)

        if np.isclose(delta_1, delta_2):
            # Hm. the distances are not changing
            self.delta = 0.0
            logger.debug("end set_delta: the distances are not changing")
            return

        if delta_1 < delta_2:
            self.start = middle_angle(sorted_areas[1][1].max_angle,
                                      sorted_areas[2][1].min_angle)
        else:
            self.start = middle_angle(sorted_areas[0][1].max_angle,
                                      sorted_areas[1][1].min_angle)
        logger.debug("end set_delta: delta wechselt: delta={}"
                     .format(self.delta))

    def symmetry_lines(self, startangle, endangle):
        logger.debug("begin symmetry_lines of {} ({}, {})"
                     .format(self.id,
                             startangle,
                             endangle))
        if less_equal(endangle, startangle):
            endangle += 2*np.pi

        angle = self.start
        while less(angle, startangle):
            angle += self.delta
        while greater(angle, startangle+self.delta):
            angle -= self.delta

        # Damit man anschliessend ohne Umstände schneiden kann.
        self.sym_startangle = angle
        self.sym_endangle = angle + self.delta
        logger.debug(" - delta: {}, sym start: {}, end: {}"
                     .format(self.sym_startangle,
                             self.delta,
                             self.sym_endangle))
        while angle < endangle:
            yield angle
            angle += self.delta
        logger.debug("end symmetry_lines")

    def minmax(self):
        mm = [99999, -99999, 99999, -99999]

        for e in self.area:
            n = e.minmax()
            mm[0] = min(mm[0], n[0])
            mm[1] = max(mm[1], n[1])
            mm[2] = min(mm[2], n[2])
            mm[3] = max(mm[3], n[3])
        return mm

    def intersect_line(self, line):
        for e in self.area:
            if e.intersect_line(line):
                return True
        return False

    def get_point_inside(self, geom):
        """return point inside area"""
        mm = self.minmax()
        y = (mm[2]+mm[3])/2
        p1 = (mm[0]-5, y)
        p2 = (mm[1]+5, y)
        line = Line(Element(start=p1, end=p2))

        points = []
        for e in self.area:
            points += e.intersect_line(line, geom.rtol, geom.atol, True)

        if len(points) < 2:
            logger.debug("WARNING: get_point_inside() failed ({})".
                         format(len(points)))
            return None

        assert(len(points) > 1)

        my_points_sorted = [(p[0], p) for p in points]
        my_points_sorted.sort()
        my_p1 = my_points_sorted[0][1]   # Startpoint
        my_p2 = my_points_sorted[-1][1]  # Endpoint

        all_points_sorted = []
        for e in geom.elements(Shape):
            points = e.intersect_line(line, geom.rtol, geom.atol, True)
            for p in points:
                if greater(p[0], my_p1[0], rtol=1e-8):
                    if less(p[0], my_p2[0], rtol=1e-8):
                        all_points_sorted.append((p[0], p))

        if len(all_points_sorted) == 0:
            p_inside = ((my_p1[0]+my_p2[0])/2, y)
            return p_inside

        all_points_sorted.sort()
        all_p1 = all_points_sorted[0][1]
        all_p2 = all_points_sorted[-1][1]
        d1 = all_p1[0] - my_p1[0]
        d2 = my_p2[0] - all_p2[0]
        if d1 > d2:
            p_inside = ((my_p1[0]+all_p1[0])/2, y)
        else:
            p_inside = ((my_p2[0]+all_p2[0])/2, y)
        return p_inside

    def render(self, renderer, color='black', with_nodes=False, fill=True):
        if fill:
            if self.render_fill(renderer):
                color = 'black'

        for e in self.area:
            e.render(renderer, color, with_nodes)
        return

    def render_fill(self, renderer, alpha=1.0):
        color = self.color()
        if not color:
            return False

        if self.is_circle():
            e = self.area[0]
            renderer.fill_circle(e.center, e.radius, color, alpha)
        else:
            nodes = [n for n in self.virtual_nodes()]
            x = [n[0] for n in nodes]
            y = [n[1] for n in nodes]
            renderer.fill(x, y, color, alpha)
        return True

    def render_legend(self, renderer, alpha=1.0):
        return renderer.new_legend_handle(self.color(), alpha, self.name())

    def remove_edges(self, g, ndec):
        for e in self.area:
            try:
                g.remove_edge(e.node1(ndec), e.node2(ndec))
            except Exception:
                continue

    def is_circle(self):
        e = self.area[0]
        if len(self.area) == 1:
            return isinstance(e, Circle) and not isinstance(e, Arc)

        if isinstance(e, Arc):
            c = e.center
            r = e.radius
            a = 0.0
            for e in self.area:
                if not isinstance(e, Arc):
                    return False
                if not points_are_close(c, e.center):
                    return False
                if not np.isclose(r, e.radius):
                    return False
                a += e.get_angle_of_arc()
            return np.isclose(a, 2.0*np.pi)

        return False

    def is_half_circle(self, center, angle):
        for e in self.area:
            if isinstance(e, Line):
                if not np.isclose(angle, alpha_line(center, e.p1)):
                    return False
                if not np.isclose(angle, alpha_line(center, e.p2)):
                    return False
            elif isinstance(e, Arc):
                if not np.isclose(angle, alpha_line(center, e.center)):
                    return False
            else:
                return False
        return True

    def is_rectangle(self):
        lines = [[c, e.m(99999.0), e.length()]
                 for c, e in enumerate(self.area)
                 if isinstance(e, Line)]
        lines.sort()

        line_count = 1
        m_first = 0.0
        m_prev = 999.999999
        c_prev = -99
        m_all = []
        for c, m, l in lines:
            if c_prev >= 0:
                if np.isclose(m_prev, m, atol=0.001):
                    if c_prev+1 != c:
                        # Gleiche Steigung, aber keine Verlängerung
                        line_count += 1
                        m_all.append(m_prev)
                else:
                    line_count += 1
                    m_all.append(m_prev)
            else:
                m_first = m

            m_prev = m
            c_prev = c

        m_all.append(m_prev)

        if np.isclose(m_prev, m_first, atol=0.001):
            line_count -= 1

        if line_count == 4:
            logger.debug("is_rectangle: m={}".format(m_all))
            if not np.isclose(m_all[0], m_all[2], atol=0.001):
                return False
            if not np.isclose(m_all[1], m_all[3], atol=0.001):
                return False
            return True

        return False

    def is_mag_rectangle(self):
        lines_ceml = [[c, e, e.m(99999.0), e.length()]
                      for c, e in enumerate(self.area)]
        # if isinstance(e, Line)]
        if len(lines_ceml) < 4:
            return False

        logger.debug("=== BEGIN OF is_mag_rectangle() [{} lines]"
                     .format(len(lines_ceml)))

        c_prev = lines_ceml[0][0]
        a_prev = 999
        p = None

        e0 = lines_ceml[0][1]
        e0_p1 = e0.p1
        e0_p2 = e0.p2
        l_prev = lines_ceml[0][3]
        m_prev = lines_ceml[0][2]

        e1 = lines_ceml[1][1]
        e1_p1 = e1.p1
        e1_p2 = e1.p2

        if points_are_close(e0_p2, e1_p1) or \
           points_are_close(e0_p2, e1_p2):
            a_prev = alpha_line(e0_p1, e0_p2)
            p = e0_p2
        else:
            if points_are_close(e0_p1, e1_p1) or \
               points_are_close(e0_p1, e1_p2):
                a_prev = alpha_line(e0_p2, e0_p1)
                p = e0_p1
            else:
                logger.error("ERROR: is_mag_rectangle(): points are not close together")
                logger.error("       e0 p1={}, p2={}".format(e0_p1, e0_p2))
                logger.error("       e1 p1={}, p2={}".format(e1_p1, e1_p2))
                return False

        def alpha_current(p, e):
            if points_are_close(p, e.p1):
                return e.p2, alpha_line(e.p1, e.p2)
            if points_are_close(p, e.p2):
                return e.p1, alpha_line(e.p2, e.p1)
            logger.error("ERROR: is_mag_rectangle(): points are not close together")
            logger.error("       p={}, p1={}, p2={}".format(p, e.p1, e.p2))
            return None, None

        lines_clam = []
        for c, e, m, l in lines_ceml[1:]:
            p, a_curr = alpha_current(p, e)
            if not p:
                return False

            if is_same_angle(a_prev, a_curr, atol=0.01):
                # its the same angle
                # assert(np.isclose(m_prev, m, atol=0.001))
                if c_prev+1 != c:
                    logger.debug(" - ok, but not an extension")
                    # ..., but not an extension
                    lines_clam.append([c_prev, l_prev, a_prev, m])
                    l_prev = e.length()
                else:
                    # ... and an extension
                    l_prev += e.length()
                    logger.debug(" - ok, it's an extension")
            else:
                # it's a different angle
                logger.debug(" - diff, angle {} and {} not equal "
                             .format(a_prev, a_curr))
                lines_clam.append([c_prev, l_prev, a_prev, m_prev])
                l_prev = e.length()

            a_prev = a_curr
            m_prev = m
            c_prev = c

        lines_clam.append([c_prev, l_prev, a_prev, m_prev])
        if np.isclose(lines_clam[0][2], lines_clam[-1][2], atol=0.001):
            # Gleicher Winkel am Anfang und am Ende
            lines_clam[0][1] += lines_clam[-1][1]  # length
            del lines_clam[-1]
            logger.debug(" > last entry deleted")

        if len(lines_clam) < 4:
            logger.debug("=== END OF is_mag_rectangle(): NO RECTANGLE #1")
            return False

        lines_lmc = [[l, m, c] for c, l, a, m in lines_clam]
        lines_lmc.sort(reverse=True)

        if not np.isclose(lines_lmc[0][1], lines_lmc[1][1], atol=0.001):
            # Die Steigungen der zwei längsten Linien müssen gleich sein
            logger.debug("=== END OF is_mag_rectangle(): NO RECTANGLE #2")
            return False

        def excursion_to_same_direction(clam):
            if len(clam) < 4:
                return False

            alpha = alpha_angle(clam[0][2], clam[1][2])
            clockwise = not alpha < np.pi

            angle_prev = clam[1][2]
            for c, l, angle_curr, m in clam[2:]:
                alpha = alpha_angle(angle_prev, angle_curr)
                if clockwise:
                    if alpha < np.pi:
                        return False
                else:
                    if alpha > np.pi:
                        return False
                angle_prev = angle_curr
            return True  # end of all_lines_with_same_direction()

        lines_cm = [[c, m] for l, m, c in lines_lmc[0:4]]
        lines_cm.sort()

        if np.isclose(lines_cm[0][1], lines_cm[2][1], atol=0.001):
            ok = excursion_to_same_direction(lines_clam)
            logger.debug("=== END OF is_mag_rectangle(): OK = {} #1"
                         .format(ok))
            return ok
        if np.isclose(lines_cm[1][1], lines_cm[3][1], atol=0.001):
            ok = excursion_to_same_direction(lines_clam)
            logger.debug("=== END OF is_mag_rectangle(): OK = {} #2"
                         .format(ok))
            return ok

        logger.debug("=== END OF is_mag_rectangle(): NO RECTANGLE #3")
        return False

    def get_mag_orient_rectangle(self):
        lines = [[e.m(99999.0), e.length(), alpha_line(e.p1, e.p2)]
                 for e in self.area
                 if isinstance(e, Line)]
        lines.sort()

        m_prev = 999.999999
        a_prev = 0.0
        l_total = 0.0
        line_length = []
        for m, l, a in lines:
            if np.isclose(m_prev, m):
                l_total += l
            else:
                if l_total > 0.0:
                    line_length.append((l_total, m_prev, a_prev))
                l_total = l
                m_prev = m
                a_prev = a
        line_length.sort(reverse=True)

        alpha = line_length[0][2]
        if alpha < 0.0:
            alpha += np.pi
        return alpha + np.pi/2

    def get_mag_orientation(self):
        if self.mag_rectangle:
            return self.get_mag_orient_rectangle()

        if self.close_to_endangle:
            if self.close_to_startangle:
                return middle_angle(self.min_angle, self.max_angle)
            else:
                return self.max_angle
        else:
            return middle_angle(self.min_angle, self.max_angle)

    def around_windings(self, areas):
        for a in areas:
            if a.is_winding():
                if not self.is_identical(a):
                    if self.is_inside(a):
                        return True
                    elif self.is_touching(a):
                        return True
        return False

    def mark_stator_subregions(self, is_inner, mirrored, alpha,
                               center, r_in, r_out):
        alpha = round(alpha, 6)

        if self.is_circle():
            self.type = 0  # air
            return self.type

        if is_inner:
            self.close_to_ag = np.isclose(r_out, self.max_dist)
            close_to_opposition = np.isclose(r_in, self.min_dist)
            airgap_radius = r_out
            opposite_radius = r_in
            airgap_toleranz = -(self.max_dist - self.min_dist) / 50.0  # 2%
        else:
            self.close_to_ag = np.isclose(r_in, self.min_dist)
            close_to_opposition = np.isclose(r_out, self.max_dist)
            airgap_radius = r_in
            opposite_radius = r_out
            airgap_toleranz = (self.max_dist - self.min_dist) / 50.0  # 2%

        self.close_to_startangle = np.isclose(self.min_angle, 0.0,
                                              1e-04, 1e-04)
        self.close_to_endangle = np.isclose(self.max_angle, alpha,
                                            1e-04, 1e-04)

        logger.debug("\n***** mark_stator_subregions [{}] *****"
                     .format(self.id))
        logger.debug(" - close_to_ag        : %s", self.close_to_ag)
        logger.debug(" - close_to_opposition: %s", close_to_opposition)
        logger.debug(" - airgap_radius      : %3.12f", airgap_radius)
        logger.debug(" - airgap_toleranz    : %3.12f", airgap_toleranz)
        logger.debug(" - opposite radius    : %3.12f", opposite_radius)
        logger.debug(" - close_to_startangle: %s", self.close_to_startangle)
        logger.debug(" - close_to_endangle  : %s", self.close_to_endangle)
        logger.debug(" - alpha              : %3.12f", alpha)
        logger.debug(" - min_angle          : %3.12f", self.min_angle)
        logger.debug(" - max_angle          : %3.12f", self.max_angle)
        logger.debug(" - min_dist           : %3.12f", self.min_dist)
        logger.debug(" - max_dist           : %3.12f", self.max_dist)

        if close_to_opposition:
            self.type = 5  # iron yoke (Joch)
            logger.debug("***** iron yoke #1\n")
            return self.type

        if self.close_to_startangle and self.close_to_endangle:
            self.type = 5  # iron yoke (Joch)
            logger.debug("***** iron yoke #2\n")
            return self.type

        if self.close_to_ag:  # close to airgap
            mm = self.minmax_angle_dist_from_center(center,
                                                    airgap_radius +
                                                    airgap_toleranz)
            self.min_air_angle = mm[0]
            self.max_air_angle = mm[1]
            air_alpha = round(alpha_angle(mm[0], mm[1]), 3)
            logger.debug(" - min_air_alpha      : {}".format(mm[0]))
            logger.debug(" - max_air_alpha      : {}".format(mm[1]))
            logger.debug(" - air_alpha          : {}".format(air_alpha))

            if self.alpha / air_alpha > 2:
                logger.debug("***** windings near airgap\n")
                self.type = 2  # windings
            else:
                self.type = 9  # air or iron near windings and near airgap?
                logger.debug("***** air or iron ??\n")
            return self.type

        if self.close_to_startangle:
            if self.is_half_circle(center, self.min_angle):
                self.type = 0  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        if self.close_to_endangle:
            if self.is_half_circle(center, self.max_angle):
                self.type = 0  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        if self.min_angle > 0.001:
            if self.max_angle < alpha - 0.001:
                self.type = 2  # windings
                logger.debug("***** windings #1\n")
            elif mirrored:
                self.type = 2  # windings
                logger.debug("***** windings #2\n")
            else:
                self.type = 0  # air
                logger.debug("***** air #2\n")
            return self.type

        logger.debug("***** air #3\n")
        return 0

    def mark_rotor_subregions(self, is_inner, mirrored, alpha,
                              center, r_in, r_out):
        logger.debug("mark_rotor_subregions")

        alpha = round(alpha, 6)

        if self.is_circle():
            self.type = 0  # air
            logger.debug(">>> air is a circle")
            return self.type

        if is_inner:
            self.close_to_ag = np.isclose(r_out, self.max_dist, atol=0.005)
            close_to_opposition = greater_equal(r_in * 1.05, self.min_dist)
            airgap_radius = r_out
            opposite_radius = r_in
            airgap_toleranz = -(self.max_dist - self.min_dist) / 50.0  # 2%
        else:
            self.close_to_ag = np.isclose(r_in, self.min_dist, atol=0.005)
            close_to_opposition = greater_equal(self.max_dist * 1.05, r_out)
            airgap_radius = r_in
            opposite_radius = r_out
            airgap_toleranz = (self.max_dist - self.min_dist) / 50.0  # 2%

        self.close_to_startangle = np.isclose(self.min_angle, 0.0)
        self.close_to_endangle = np.isclose(self.max_angle, alpha)

        logger.debug("\n***** mark_rotor_subregions [{}] *****"
                     .format(self.id))
        logger.debug(" - close_to_ag        : %s", self.close_to_ag)
        logger.debug(" - close_to_opposition: %s", close_to_opposition)
        logger.debug(" - min dist           : %3.12f", self.min_dist)
        logger.debug(" - max dist           : %3.12f", self.max_dist)
        logger.debug(" - airgap radius      : %3.12f", airgap_radius)
        logger.debug(" - opposite radius    : %3.12f", opposite_radius)
        logger.debug(" - close_to_startangle: %s", self.close_to_startangle)
        logger.debug(" - close_to_endangle  : %s", self.close_to_endangle)
        logger.debug(" - alpha              : %3.12f", alpha)
        logger.debug(" - min_angle          : %3.12f", self.min_angle)
        logger.debug(" - max_angle          : %3.12f", self.max_angle)

        if close_to_opposition:
            self.type = 1  # iron
            logger.debug("***** iron (close to opposition)\n")
            return self.type

        self.mag_rectangle = self.is_mag_rectangle()

        if self.close_to_ag:
            mm = self.minmax_angle_dist_from_center(center,
                                                    airgap_radius +
                                                    airgap_toleranz)
            air_alpha = round(alpha_angle(mm[0], mm[1]), 3)
            logger.debug(" - air_alpha          : {}".format(air_alpha))

            if air_alpha / alpha < 0.2:
                self.phi = self.get_mag_orientation()
                self.type = 8  # air or magnet ?
                logger.debug("***** air #1 (close to airgap)\n")
                return self.type

            if air_alpha / alpha > 0.6:
                self.phi = self.get_mag_orientation()
                self.type = 3  # magnet
                logger.debug("***** magnet (close to airgap)\n")
            else:
                self.phi = self.get_mag_orientation()
                self.type = 9  # iron or magnet ?
                logger.debug("***** iron or magnet(close to airgap)\n")
            return self.type

        if self.is_mag_rectangle():
            self.phi = self.get_mag_orientation()
            self.type = 4  # magnet embedded
            logger.debug("***** magnet (embedded)\n")
            return self.type

        if not (self.close_to_startangle or self.close_to_endangle):
            self.type = 0  # air
            logger.debug("***** air (somewhere)\n")
            return self.type

        if self.close_to_startangle:
            if self.is_half_circle(center, self.min_angle):
                self.type = 0  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        if self.close_to_endangle:
            if self.is_half_circle(center, self.max_angle):
                self.type = 0  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        self.type = 0  # air
        logger.debug("***** air (remains)\n")
        return self.type

    def mark_unknown_subregions(self, mirrored, alpha,
                                center, r_in, r_out):
        logger.debug("mark_unknown_subregions")

        if self.is_circle():
            self.type = 0  # air
            logger.debug(">>> air is a circle")
            return self.type

        if self.is_mag_rectangle():
            self.type = 4  # magnet embedded
            logger.debug(">>> magnet embedded")
            self.phi = self.get_mag_orient_rectangle()
            return self.type

        close_to_max_radius = np.isclose(r_out, self.max_dist)
        close_to_min_radius = np.isclose(r_in, self.min_dist)

        if close_to_max_radius and close_to_min_radius:
            self.type = 1  # iron
            logger.debug(">>> iron close to min- and max-radius")
            return self.type

        self.close_to_startangle = np.isclose(self.min_angle, 0.0)
        self.close_to_endangle = np.isclose(self.max_angle, alpha)

        if self.close_to_startangle and self.close_to_endangle:
            self.type = 1  # iron
            logger.debug(">>> iron close to start- and end-angle")
            return self.type

        self.type = 0  # air
        logger.debug(">>> air remains")
        return self.type

    def print_area(self):
        center = [0.0, 0.0]
        for s in self.area:
            mm = s.minmax_angle_from_center(center)
            print(" --- angle min={}, max={}".format(mm[0], mm[1]))

    def __lt__(self, a):
        if self.symmetry != a.symmetry:
            return self.symmetry > a.symmetry

        if self.sym_type != a.sym_type:
            return self.sym_type > a.sym_type

        if self.count != a.count:
            return self.count > a.count

        if not np.isclose(self.height, a.height):
            return self.height > a.height

        if self.sym_tolerance > 0.0:
            if not np.isclose(round(self.min_dist, 4),
                              round(a.min_dist, 4), 1e-03,
                              self.sym_tolerance):
                return less_equal(self.min_dist, a.min_dist)
            if not np.isclose(round(self.max_dist, 4),
                              round(a.max_dist, 4), 1e-03,
                              self.sym_tolerance):
                return less_equal(self.max_dist, a.max_dist)
            if not np.isclose(round(self.alpha, 2),
                              round(a.alpha, 2), 1e-01, 1e-01):
                return less_equal(self.alpha, a.alpha)
        else:
            if not np.isclose(round(self.min_dist, 2),
                              round(a.min_dist, 2)):
                return less_equal(self.min_dist, a.min_dist)
            if not np.isclose(round(self.max_dist, 2),
                              round(a.max_dist, 2)):
                return less_equal(self.max_dist, a.max_dist)
            if not np.isclose(round(self.alpha, 2),
                              round(a.alpha, 2), 1e-01, 1e-02):
                return less_equal(self.alpha, a.alpha)

        return self.min_angle < a.min_angle

    def __str__(self):
        return "Area {}\n".format(self.id) + \
            "distance: from {} to {}\n".\
            format(round(self.min_dist, 4), round(self.max_dist, 4)) + \
            "height..: {}\n".format(self.height) + \
            "alpha...: {}\n".format(self.alpha) + \
            "angle...: from {} to {}\n".\
            format(round(self.min_angle, 6), round(self.max_angle, 6)) + \
            "delta...: {}\n".format(self.delta) + \
            "number..: {}\n".format(self.count) + \
            "equal...: {}\n".format(len(self.equal_areas)) + \
            "symmetry: {}\n".format(self.symmetry) + \
            "sym_type: {}".format(self.sym_type)
