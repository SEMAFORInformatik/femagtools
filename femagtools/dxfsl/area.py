# -*- coding: utf-8 -*-
"""manage a geometry with
    lines, circles and arcs built from DXF

  NOTE: This code is in highly experimental state.
        Use at your own risk.

  Author: Ronald Tanner
    Date: 2017/07/06
"""
from __future__ import print_function
import numpy as np
import networkx as nx
import logging
from .functions import less_equal, less, greater_equal, greater
from .functions import distance, alpha_angle, alpha_line, min_angle, max_angle
from .functions import point, line_m, line_n, intersect_point, points_are_close
from .functions import middle_angle, part_of_circle
from .shape import Element, Shape, Line, Arc, Circle

logger = logging.getLogger('femagtools.area')


#############################
#            Area           #
#############################

class Area(object):
    def __init__(self, area, center, sym_tolerance):
        self.area = area
        self.type = 0  # material
        self.phi = 0.0
        self.min_angle = 0.0
        self.max_angle = 0.0
        self.close_to_startangle = False
        self.close_to_endangle = False
        self.min_dist = 99999.0
        self.max_dist = 0.0
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
            return 'blue'
        if self.type == 6:
            return 'darkblue'
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

            mm_angle = s.minmax_angle_from_center(center)
            self.min_angle = min_angle(self.min_angle, mm_angle[0])
            self.max_angle = max_angle(self.max_angle, mm_angle[1])

        self.alpha = round(alpha_angle(self.min_angle, self.max_angle), 3)

    def minmax_angle_dist_from_center(self, center, dist):
        s = self.area[0]
        my_min_angle = self.max_angle
        my_max_angle = self.min_angle
        mm_angle = None
        for s in self.area:
            mm_angle = s.minmax_angle_dist_from_center(center, dist)
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

    def get_most_left_point(self, center, radius, angle):
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

        return (dist, the_axis_p, the_area_p)

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
                print(" - OK")
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
        self.delta = 0.0
        self.symmetry = 0

        if len(self.equal_areas) < 2:
            # Mit zwei Objekten lässt sich das Teil nur noch halbieren. Das
            # wird zum Schluss sowieso versucht.
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

        if len(delta_sorted) == 1:
            # simple case: all have the same angle
            self.delta = alpha_angle(sorted_areas[0][1].min_angle,
                                     sorted_areas[1][1].min_angle)
            self.start = middle_angle(sorted_areas[0][1].max_angle,
                                      sorted_areas[1][1].min_angle)
            self.sym_type = 3
            self.symmetry = part_of_circle(0.0, self.delta, 1)
            return

        if len(delta_sorted) > 2:
            # Mehr als 2 Winkel untersuchen wir (noch) nicht. Wir brechen
            # die Suche nach dem richtigen Winkel ab.
            return

        # Bei 2 verschiedenen Winkeln werden die näher beieinander liegenden
        # Objekte zusammen genommen.

        if len(self.equal_areas) < 4:
            # Wenn nicht mehr als 4 Objekte vorhanden sind, brechen wir auch
            # ab.
            return

        percent = delta_sorted[0][0] / (len(self.equal_areas)+1)
        if percent > 0.75:
            # lets assume we only have on angle
            self.delta = alpha_angle(sorted_areas[0][1].min_angle,
                                     sorted_areas[1][1].min_angle)
            self.start = middle_angle(sorted_areas[0][1].max_angle,
                                      sorted_areas[1][1].min_angle)
            self.sym_type = 2
            self.symmetry = part_of_circle(0.0, self.delta, 1)
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
            return

        if delta_1 < delta_2:
            self.start = middle_angle(sorted_areas[1][1].max_angle,
                                      sorted_areas[2][1].min_angle)
        else:
            self.start = middle_angle(sorted_areas[0][1].max_angle,
                                      sorted_areas[1][1].min_angle)

    def symmetry_lines(self, startangle, endangle):
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

        while angle < endangle:
            yield angle
            angle += self.delta

    def minmax(self):
        mm = [99999, -99999, 99999, -99999]

        for e in self.area:
            n = e.minmax()
            mm[0] = min(mm[0], n[0])
            mm[1] = max(mm[1], n[1])
            mm[2] = min(mm[2], n[2])
            mm[3] = max(mm[3], n[3])
        return mm

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

    def render(self, renderer, color='black', with_nodes=False):
        for e in self.area:
            e.render(renderer, color, with_nodes)
        return

    def render_fill(self, renderer, alpha=1.0):
        color = self.color()
        if not color:
            return

        if self.is_circle():
            e = self.area[0]
            renderer.fill_circle(e.center, e.radius, color, alpha)
        else:
            nodes = [n for n in self.virtual_nodes()]
            x = [n[0] for n in nodes]
            y = [n[1] for n in nodes]
            renderer.fill(x, y, color, alpha)

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

    def is_rectangle(self):
        lines = [[c, e.m(99999.0), e.length()]
                 for c, e in enumerate(self.area)
                 if isinstance(e, Line)]
        lines.sort()

        line_count = 1
        m_first = 0.0
        m_prev = 999.999999
        c_prev = -99
        for c, m, l in lines:
            if c_prev >= 0:
                if np.isclose(m_prev, m, atol=0.001):
                    if c_prev+1 != c:
                        # Gleiche Steigung, aber keine Verlängerung
                        line_count += 1
                else:
                    line_count += 1
            else:
                m_first = m
            m_prev = m
            c_prev = c

        if np.isclose(m_prev, m_first, atol=0.001):
            line_count -= 1

        return line_count == 4

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

    def around_windings(self, areas):
        for a in areas:
            if a.is_winding():
                if not self.is_identical(a):
                    if self.is_inside(a):
                        return True
        return False

    def mark_stator_subregions(self, is_inner, mirrored, alpha,
                               center, r_in, r_out):
        alpha = round(alpha, 6)

        if self.is_circle():
            self.type = 0  # air
            return self.type

        if is_inner:
            close_to_ag = np.isclose(r_out, self.max_dist)
            close_to_opposition = np.isclose(r_in, self.min_dist)
            airgap_radius = r_out
        else:
            close_to_ag = np.isclose(r_in, self.min_dist)
            close_to_opposition = np.isclose(r_out, self.max_dist)
            airgap_radius = r_in

        self.close_to_startangle = np.isclose(self.min_angle, 0.0)
        self.close_to_endangle = np.isclose(self.max_angle, alpha)

        if close_to_opposition:
            self.type = 5  # iron yoke (Joch)
            return self.type

        if self.close_to_startangle and self.close_to_endangle:
            self.type = 5  # iron yoke (Joch)
            return self.type

        if close_to_ag:  # close to airgap
            mm = self.minmax_angle_dist_from_center(center, airgap_radius)
            air_alpha = round(alpha_angle(mm[0], mm[1]), 3)

            if air_alpha / alpha < 0.2:
                self.type = 0  # air
                return self.type

            if air_alpha / alpha < 0.5:
                self.type = 9  # air or iron near windings?
            else:
                self.type = 6  # iron shaft (Zahn)
            return self.type

        if self.min_angle > 0.001:
            if self.max_angle < alpha - 0.001:
                self.type = 2  # windings
            elif mirrored:
                self.type = 2  # windings
            else:
                self.type = 0  # air
            return self.type

        return 0

    def mark_rotor_subregions(self, is_inner, mirrored, alpha,
                              center, r_in, r_out):
        my_alpha = round(self.max_angle - self.min_angle, 6)
        alpha = round(alpha, 6)

        if self.is_circle():
            self.type = 0  # air
            return self.type

        if is_inner:
            close_to_ag = np.isclose(r_out, self.max_dist)
            close_to_opposition = np.isclose(r_in, self.min_dist)
            airgap_radius = r_out
        else:
            close_to_ag = np.isclose(r_in, self.min_dist)
            close_to_opposition = np.isclose(r_out, self.max_dist)
            airgap_radius = r_in

        self.close_to_startangle = np.isclose(self.min_angle, 0.0)
        self.close_to_endangle = np.isclose(self.max_angle, alpha)

        if close_to_opposition:
            self.type = 1  # iron
            return self.type

        if close_to_ag:
            mm = self.minmax_angle_dist_from_center(center, airgap_radius)
            air_alpha = round(alpha_angle(mm[0], mm[1]), 3)
            if air_alpha / alpha < 0.2:
                self.type = 0  # air
                return self.type

            if air_alpha / alpha > 0.6:
                self.type = 3  # magnet
                if self.close_to_endangle:
                    if self.close_to_startangle:
                        self.phi = middle_angle(self.min_angle, self.max_angle)
                    else:
                        self.phi = self.max_angle
                else:
                    self.phi = middle_angle(self.min_angle, self.max_angle)
            else:
                self.type = 1  # iron
            return self.type

        if my_alpha / alpha > 0.5:
            if self.is_rectangle():
                self.type = 4  # magnet embedded
                self.phi = self.get_mag_orient_rectangle()
                return self.type

            if not (self.close_to_startangle or self.close_to_endangle):
                self.type = 0  # air
                return self.type

        self.type = 1  # iron
        if self.min_angle > 0.001:
            if my_alpha / alpha < 0.4:
                if self.max_angle < alpha - 0.001:
                    self.type = 0  # air
                elif mirrored:
                    self.type = 0  # air

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
        return "Area\ndistance: from {} to {}\n".\
            format(round(self.min_dist, 4), round(self.max_dist, 4)) + \
            "alpha...: {}\n".format(self.alpha) + \
            "angle...: from {} to {}\n".\
            format(round(self.min_angle, 6), round(self.max_angle, 6)) + \
            "delta...: {}".format(self.delta)
