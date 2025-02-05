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
from .functions import greater_angle, less_angle
from .functions import distance, alpha_angle, alpha_line, min_angle, max_angle
from .functions import point, line_m, line_n, intersect_point, points_are_close
from .functions import middle_angle, part_of_circle, is_same_angle
from .functions import area_size, normalise_angle, positive_angle
from .shape import Element, Shape, Line, Arc, Circle, is_Circle, is_Line, is_Arc

logger = logging.getLogger('femagtools.area')


#############################
#            Area           #
#############################

area_number = 0

TYPE_AIR = 0
TYPE_IRON = 1
TYPE_WINDINGS = 2
TYPE_MAGNET_AIRGAP = 3
TYPE_MAGNET_RECT = 4
TYPE_YOKE = 5
TYPE_TOOTH = 6
TYPE_MAGNET_OR_AIR = 8
TYPE_AIR_OR_IRON = 9
TYPE_MAGNET_OR_IRON = 9
TYPE_SHAFT = 10
TYPE_MAGNET_RECT_NEAR_AIRGAP = 11
TYPE_WINDINGS_OR_AIR = 12
TYPE_WINDINGS_OR_IRON = 13
TYPE_FD_WINDINGS = 15
TYPE_MAGNET_UNDEFINED = 99
TYPE_GROUP = 20

class Area(object):
    def __init__(self, area, center, sym_tolerance):
        self.area = area
        self.type = -1  # material
        self.phi = 0.0
        self.min_angle = 0.0
        self.max_angle = 0.0
        self.min_air_angle = 0.0
        self.max_air_angle = 0.0
        self.close_to_ag = False
        self.close_to_ag_startcorner = False
        self.close_to_ag_endcorner = False
        self.close_to_startangle = False
        self.close_to_endangle = False
        self.mag_rectangle = False
        self.mag_width = 0.0
        self.min_dist = 99999.0
        self.max_dist = 0.0
        self.min_x = None
        self.max_x = None
        self.min_y = None
        self.max_y = None
        self.height = 0.0
        self.alpha = 0.0
        self.count = 1
        self.equal_areas = []
        self.delta = 0.0
        self.start = 0.0
        self.sym_startangle = 0.0
        self.sym_endangle = 0.0
        self.sym_upper_left_dist = None
        self.sym_upper_right_dist = None
        self.sym_lower_left_dist = None
        self.sym_lower_right_dist = None
        self.sym_type = 0
        self.symmetry = 0
        self.sym_tolerance = sym_tolerance
        self.calc_signature(center)
        self.surface = 0.0
        global area_number
        area_number += 1
        self.id = area_number
        self.is_child = False
        self.areas_inside = {}
        self.areas_of_group = []
        self.group_is_inside = False

    def identifier(self):
        return "{}-{}".format(self.id, self.type)

    def get_id(self):
        return self.id

    def number_of_elements(self):
        return len(self.area)

    def elements(self):
        return self.area

    def copy_of_elements(self):
        return [e.clone() for e in self.elements() if e]

    def list_of_nodes(self):
        if len(self.area) < 1:
            return

        e0 = self.area[0]
        if len(self.area) == 1:
            yield e0.n1
            yield e0.n2
            return

        e1 = self.area[1]
        try:
            if e1.get_node_number(e0.n1, override=True) == 0:
                nx = e0.n2
            else:
                nx = e0.n1
            yield nx

            for e1 in self.area[1:]:
                if e1.get_node_number(nx) == 1:
                    nx = e1.n2
                else:
                    nx = e1.n1
                yield nx
        except ValueError as e:
            logger.error("list_of_nodes(): FATAL ERROR: %s", e)
            return
        except Exception:
            return

    def list_of_elements(self):
        if len(self.area) < 2:
            return

        e0 = self.area[0]
        e1 = self.area[1]
        try:
            if e1.get_node_number(e0.n1, override=True) == 0:
                n1 = e0.n1
                n2 = e0.n2
            else:
                n1 = e0.n2
                n2 = e0.n1
            yield n1, n2, e0

            for e1 in self.area[1:]:
                if e1.get_node_number(n2) == 1:
                    n1 = e1.n1
                    n2 = e1.n2
                else:
                    n1 = e1.n2
                    n2 = e1.n1
                yield n1, n2, e1
        except ValueError as e:
            logger.error("list_of_elements(): FATAL ERROR: %s", e)
            return
        except Exception:
            return

    def has_fsl(self):
        try:
            return self.fsl
        except AttributeError:
            return True

    def list_of_equal_edges(self, a):
        for e1 in self.area:
            for e2 in a.area:
                if e1.n1 == e2.n1 and e1.n2 == e2.n2:
                    yield e1

    def reduce_line_nodes(self, geom, mindist=0.01):
        """reduces number of nodes (lines only)
        https://rdp.readthedocs.io/en/latest
        Note: this feature is deactivated silently
          if the rdp package is not installed.
        """
        try:
            import rdp
        except ModuleNotFoundError:
            return 0

        def is_valid_line_(n1, e):
            if not isinstance(e, Line):
                return False
            if e.has_attribute('del'):
                return False
            return True

        def reduce_nodes_(lines, mindist):
            if not len(lines) > 1:
                return 0

            nodes = [n1 for n1, n2, e in lines]
            n1, n2, e = lines[-1]
            nodes.append(n2)

            remaining_nodes = rdp.rdp(nodes,
                                      epsilon=mindist)
            nodes_deleted = len(nodes) - len(remaining_nodes)
            if not nodes_deleted:
                return 0

            for n1, n2, e in lines:
                e.set_attribute('del')
                e.set_my_color('yellow')
                geom.remove_edge(e)

            n1 = remaining_nodes[0]
            for n2 in remaining_nodes[1:]:
                geom.add_line(n1, n2)
                n1 = n2

            self.area = []
            return nodes_deleted

        # -----
        nodes_deleted = 0
        lines = []
        for n1, n2, e in self.list_of_elements():
            if not is_valid_line_(n1, e):
                nodes_deleted += reduce_nodes_(lines, mindist)
                lines = []
            elif not geom.num_of_neighbors(n1) == 2:
                nodes_deleted += reduce_nodes_(lines, mindist)
                lines = [(n1, n2, e)]
            else:
                lines.append((n1, n2, e))

        nodes_deleted += reduce_nodes_(lines, mindist)
        return nodes_deleted

    def reduce_element_nodes(self, geom, mindist=0.01):
        """reduces number of nodes (lines only)
        https://rdp.readthedocs.io/en/latest
        Note: this feature is deactivated silently
          if the rdp package is not installed.
        """
        corners = geom.start_corners + geom.end_corners

        try:
            import rdp
        except ModuleNotFoundError:
            return 0

        def is_valid_element_(n1, e):
            return not e.has_attribute('del')

        def reduce_nodes_(elements, mindist):
            if not len(elements) > 1:
                return 0

            old_nodes = []
            for n1, n2, e in elements:
                old_nodes.append(n1)
                nodes = [n for n in e.get_nodes(parts=24)]
                if len(nodes) > 2:
                    if points_are_close(n1, nodes[0]):
                        old_nodes += nodes[1:-1]
                    elif points_are_close(n2, nodes[0]):
                        nodes.reverse()
                        old_nodes += nodes[1:-1]
            n1, n2, e = elmts[-1]
            old_nodes.append(n2)

            new_nodes = rdp.rdp(old_nodes,
                                epsilon=mindist)
            nodes_deleted = len(old_nodes) - len(new_nodes)
            if not nodes_deleted:
                return 0

            for n1, n2, e in elements:
                e.set_attribute('del')
                e.set_my_color('yellow')
                geom.remove_edge(e)

            n1 = new_nodes[0]
            for n2 in new_nodes[1:]:
                geom.add_line(n1, n2)
                n1 = n2

            self.area = []
            return nodes_deleted

        # -----
        nodes_deleted = 0
        tiny_mindist = 0.05
        elmts = []
        has_tiny = False
        for n1, n2, e in self.list_of_elements():
            if not is_valid_element_(n1, e):
                if has_tiny:
                    nodes_deleted += reduce_nodes_(elmts, mindist)
                elmts = []
            elif not geom.num_of_neighbors(n1) == 2 or n1 in corners:
                if has_tiny:
                    nodes_deleted += reduce_nodes_(elmts, mindist)
                has_tiny = e.is_tiny(tiny_mindist)
                elmts = [(n1, n2, e)]
            else:
                if e.is_tiny(tiny_mindist):
                    has_tiny = True
                elmts.append((n1, n2, e))
        if has_tiny:
            nodes_deleted += reduce_nodes_(elmts, mindist)
        return nodes_deleted

    def virtual_nodes(self, render=False, parts=64):
        if len(self.area) < 2:
            return

        prev_nodes = [n for n in self.area[0].get_nodes(parts=parts,
                                                        render=render)]
        next_nodes = [n for n in self.area[1].get_nodes(parts=parts,
                                                        render=render)]
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
            next_nodes = [n for n in e.get_nodes(parts=parts, render=render)]

            if points_are_close(next_nodes[-1], last_point, 1e-03, 1e-01):
                next_nodes = next_nodes[::-1]
            for n in next_nodes:
                yield n
            last_point = next_nodes[-1]

    def legend(self):
        if self.type == TYPE_IRON:
            return 'Iron'
        if self.type == TYPE_WINDINGS:
            return 'Windings'
        if self.type == TYPE_FD_WINDINGS:
            return 'Field Windings'
        if self.type == TYPE_MAGNET_AIRGAP or self.type == TYPE_MAGNET_RECT:
            return 'Magnet'
        if self.type == TYPE_YOKE:
            return 'Yoke'
        if self.type == TYPE_TOOTH:
            return 'Tooth'
        if self.type == TYPE_SHAFT:
            return 'Shaft'
        return ''

    def name(self):
        if self.type == TYPE_IRON:
            return 'Iron'
        if self.type == TYPE_WINDINGS:
            return 'Wndg'
        if self.type == TYPE_FD_WINDINGS:
            return 'FD_Wndg'
        if self.type == TYPE_MAGNET_AIRGAP or self.type == TYPE_MAGNET_RECT:
            return 'Mag'
        if self.type == TYPE_YOKE:
            return 'StJo'
        if self.type == TYPE_TOOTH:
            return 'StZa'
        if self.type == TYPE_SHAFT:
            return 'Shft'
        return ''

    def color(self):
        if self.type == TYPE_IRON:
            return 'cyan'
        if self.type == TYPE_WINDINGS:
            return 'green'
        if self.type == TYPE_FD_WINDINGS:
            return 'yellow'
        if self.type == TYPE_MAGNET_AIRGAP or self.type == TYPE_MAGNET_RECT:
            return 'red'
        if self.type == TYPE_YOKE:
            return 'cyan'
        if self.type == TYPE_TOOTH:
            return 'skyblue'
        if self.type == TYPE_SHAFT:
            return 'lightgrey'
        return 'white'

    def color_alpha(self):
        if self.type == TYPE_IRON:
            return 0.3
        if self.type == TYPE_WINDINGS:
            return 1.0
        if self.type == TYPE_FD_WINDINGS:
            return 1.0
        if self.type == TYPE_MAGNET_AIRGAP or self.type == TYPE_MAGNET_RECT:
            return 1.0
        if self.type == TYPE_YOKE:
            return 0.5
        if self.type == TYPE_TOOTH:
            return 1.0
        if self.type == TYPE_SHAFT:
            return 0.8
        return 1.0

    def is_iron(self):
        return \
            self.type == TYPE_IRON or \
            self.type == TYPE_YOKE or \
            self.type == TYPE_TOOTH

    def is_stator_iron_yoke(self):
        return self.type == TYPE_YOKE

    def is_stator_iron_tooth(self):
        return self.type == TYPE_TOOTH

    def is_rotor_iron(self):
        return self.type == TYPE_IRON

    def is_winding(self):
        return \
            self.type == TYPE_WINDINGS or \
            self.type == TYPE_FD_WINDINGS

    def is_field_winding(self):
        return self.type == TYPE_FD_WINDINGS

    def is_magnet(self):
        return self.type == TYPE_MAGNET_AIRGAP or self.type == TYPE_MAGNET_RECT

    def is_shaft(self):
        return self.type == TYPE_SHAFT

    def is_air(self):
        return self.type == TYPE_AIR

    def is_type(self, type):
        return self.type == type

    def set_type(self, t):
        self.type = t

    def calc_signature(self, center):
        if not self.area:
            return

        self.min_x, self.max_x, self.min_y, self.max_y = self.minmax()
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

    def center_is_inside(self, center):
        if less(self.min_x, center[0], rtol=1e-03, atol=1e-04) and \
           greater(self.max_x, center[0], rtol=1e-03, atol=1e-04) and \
           less(self.min_y, center[1], rtol=1e-03, atol=1e-04) and \
           greater(self.max_y, center[1], rtol=1e-03, atol=1e-04):
            return True
        return False

    def minmax_dist_from_center(self, center):
        nodes = np.array([n for e in self.area for n in e.get_nodes()])
        return (np.min(np.linalg.norm(nodes-center, axis=1)),
                np.max(np.linalg.norm(nodes-center, axis=1)))

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

    def is_inside(self, area, geom):
        if less_equal(area.min_dist, self.min_dist, rtol=1e-8):
            return False
        if greater_equal(area.max_dist, self.max_dist, rtol=1e-8):
            return False
        if less_equal(area.min_angle, self.min_angle, rtol=1e-8):
            return False
        if greater_equal(area.max_angle, self.max_angle, rtol=1e-8):
            return False

        p1 = self.get_point_inside(geom)
        p2 = area.get_simple_point_inside(geom)

        if p1 is None or p2 is None:
            logger.debug(" -- Line from %s to %s", p1, p2)
            logger.debug(" self has %s elements", len(self.area))
            logger.debug(" area has %s elements", len(area.area))
            return False

        line = Line(Element(start=p1, end=p2))
        plist = self.get_intersect_points(line)
        points = len(plist)

        aux_lines = [e for e in self.area if e.has_attribute('auxline')]
        if aux_lines:
            logger.debug("=> is inside = True (auxlines)")
            return True

        if plist:
            if points_are_close(p1, plist[0]):
                points -= 1
            elif points_are_close(p1, plist[-1]):
                points -= 1

        if points % 2 == 0:
            return True  # p2 is inside
        return False

    def is_touching(self, area):
        for n in self.list_of_nodes():
            x = [p for p in area.list_of_nodes() if points_are_close(n, p)]
            if x:
                return True
        return False

    def is_touching_both_sides(self):
        return (self.close_to_startangle and self.close_to_endangle)

    def is_in_touch_with_area(self, geom, a):
        n1 = self.area[0].n1
        n2 = a.area[0].n2
        try:
            return nx.has_path(geom.g, n1, n2)
        except nx.NetworkXError:
            logger.warning("has_path() failed")
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
        if a.get_id() < self.get_id():
            dist_id = '{}-{}'.format(a.get_id(), self.get_id())
        else:
            dist_id = '{}-{}'.format(self.get_id(), a.get_id())

        for p1 in self.list_of_nodes():
            for p2 in a.list_of_nodes():
                d = distance(p1, p2)
                gap_list.append((d, (p1, p2), dist_id, a.get_id()))

        if rightangle is not None:
            d, p1, p2 = a.get_nearest_point(center, radius, rightangle)
            if p1:
                gap_list.append((d, (p1, p2), dist_id, a.get_id()))
        if leftangle is not None:
            d, p1, p2 = a.get_nearest_point(center, radius, leftangle)
            if p1:
                gap_list.append((d, (p1, p2), dist_id, a.get_id()))
        if not gap_list:
            return []
        gap_list.sort()
        return [gap_list[0]]

    def get_nearest_point(self, center, radius, angle):
        axis_p = point(center, radius, angle)
        axis_m = line_m(center, axis_p)
        axis_n = line_n(center, axis_m)
        logger.debug("===== get_nearest_point in %s =====", self.identifier())

        the_area_p = None
        the_axis_p = None
        dist = 99999
        for n in self.list_of_nodes():
            p = intersect_point(n, center, axis_m, axis_n)
            d = distance(n, p)
            logger.debug("intersect point: %s", p)
            logger.debug("dist...........: %s", d)
            if d < dist:
                dist = d
                the_area_p = n
                the_axis_p = p

        logger.debug("min dist..: %s", dist)
        logger.debug("axis point: %s", the_axis_p)
        logger.debug("area point: %s", the_area_p)
        logger.debug("=============================")

        if the_area_p is None:
            return (None, None, None)

        return (dist,
                (the_axis_p[0], the_axis_p[1]),
                (the_area_p[0], the_area_p[1]))

    def get_alpha(self, center):
        if self.center_is_inside(center):
            return np.pi*2.0
        return alpha_angle(self.min_angle,
                           self.max_angle,
                           rtol=0.0,
                           atol=0.0)

    def get_mid_angle(self, center):
        if self.center_is_inside(center):
            return np.pi
        return middle_angle(self.min_angle, self.max_angle)

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

    def set_symmetry_parameter(self, center):
        all_list = [(distance(center, n), alpha_line(center, n))
                    for n in self.list_of_nodes()]
        mid = middle_angle(self.min_angle, self.max_angle)
        left_list = [(d, a) for d, a in all_list if greater_angle(a, mid)]
        right_list = [(d, a) for d, a in all_list if less_angle(a, mid)]
        left_list.sort()
        right_list.sort()

        if left_list:
            l_low_d, l_low_a = left_list[0]
            l_up_d, l_up_a = left_list[-1]
        else:
            l_low_d = self.min_dist
            l_up_d = self.max_dist
        if right_list:
            r_low_d, r_low_a = right_list[0]
            r_up_d, r_up_a = right_list[-1]
        else:
            r_low_d = self.min_dist
            r_up_d = self.max_dist
        self.sym_upper_left_dist = l_up_d
        self.sym_upper_right_dist = r_up_d
        self.sym_lower_left_dist = l_low_d
        self.sym_lower_right_dist = r_low_d

    def is_symmetry_equal(self, area):
        logger.debug("check area %s -- %s", self.get_id(), area.get_id())

        bad = False
        if not np.isclose(self.sym_lower_left_dist, area.sym_lower_left_dist,
                           rtol=5e-1, atol=5e-1):
            logger.debug("Lower left: %s != %s",
                         self.sym_lower_left_dist,
                         area.sym_lower_left_dist)
            bad = True

        if not np.isclose(self.sym_lower_right_dist, area.sym_lower_right_dist,
                          rtol=5e-1, atol=5e-1):
            logger.debug("Lower right: %s != %s",
                         self.sym_lower_right_dist,
                         area.sym_lower_right_dist)
            bad = True

        if not np.isclose(self.sym_upper_left_dist, area.sym_upper_left_dist,
                          rtol=5e-1, atol=5e-1):
            logger.debug("Upper left: %s != %s",
                         self.sym_upper_left_dist,
                         area.sym_upper_left_dist)
            bad = True

        if not np.isclose(self.sym_upper_right_dist, area.sym_upper_right_dist,
                          rtol=5e-1, atol=5e-1):
            logger.debug("Upper right: %s != %s",
                         self.sym_upper_right_dist,
                         area.sym_upper_right_dist)
            bad = True

        return not bad

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

    def intersect_area(self, line):
        for e in self.area:
            if e.intersect_shape(line, include_end=True):
                return True
        return False

    def get_intersect_points(self, line, rtol=1.e-4, atol=1e-3):
        points = []
        for e in self.area:
            points += e.intersect_line(line, rtol=rtol, atol=atol, include_end=True)
        if not points:
            return []
        points = [(x, y) for x, y in points]
        points.sort()
        p1 = points[0]
        pts = [p1]
        for p2 in points[1:]:
            if not points_are_close(p1, p2, rtol=rtol, atol=atol):
                pts.append(p2)
            p1 = p2
        return pts

    def is_point_inside(self, pt):
        for e in self.area:
            if e.is_point_inside(pt, include_end=True):
                return True
        return False

    def the_area_is_inside_area(self, a):
        p = a.area[0].n1
        return self.the_point_is_inside_area(p)

    def the_point_is_inside_area(self, p):
        x, y = p
        if less_equal(x, self.min_x) or greater_equal(x, self.max_x):
            return False
        if less_equal(y, self.min_y) or greater_equal(y, self.max_y):
            return False

        p1 = (self.min_x - 5, y)
        p2 = (self.max_x + 5, y)
        line = Line(Element(start=p1, end=p2))
        pts = self.get_intersect_points(line)
        if len(pts) % 2 != 0:
            return False

        c = 0
        for ax, ay in pts:
            if not less(ax, x):
                if c % 2 != 1:
                    return False
                break
            c += 1

        if not c < len(pts):
            return False

        ax, ay = pts[c]
        if not less(x, ax):
            return False
        return True

    def get_best_point_inside(self, geom):
        px1 = self.min_x - 5
        px2 = self.max_x + 5

        y_dist = self.max_y - self.min_y
        step = y_dist / 6
        y_list = np.arange(self.min_y + step * 0.3,
                           self.max_y - step * 0.3,
                           step)

        lines = []
        for y in y_list:
            p1 = (px1, y)
            p2 = (px2, y)
            line = Line(Element(start=p1, end=p2))
            lines.append({'line': line,
                          'pts': [],
                          'y': y,
                          'x': []})

        for e in self.area:
            points = []
            for line in lines:
                line['pts'] += e.intersect_line(line['line'],
                                                geom.rtol,
                                                geom.atol,
                                                True)
        for line in lines:
            x_sorted = [p[0] for p in line['pts']]
            x_sorted.sort()
            if x_sorted:
                line['start_x'] = x_sorted[0]
                line['end_x'] = x_sorted[-1]

        for e in geom.elements(Shape):
            for line in lines:
                if line.get('start_x', None) is None:
                    continue
                points = e.intersect_line(line['line'],
                                          geom.rtol,
                                          geom.atol,
                                          True)

                for p in points:
                    if greater(p[0],
                               line['start_x'],
                               rtol=1e-8):
                        if less(p[0],
                                line['end_x'],
                                rtol=1e-8):
                            line['x'].append(p[0])

        points = []
        for line in lines:
            if line.get('start_x', None) is None:
                continue
            line['x'].sort()
            x1 = line['start_x']
            x2 = line['end_x']
            if line['x']:
                x = line['x'][0]  # first point
                line['x_dist'] = x - x1
                points.append((line['x_dist'], (x1+x)/2, line['y']))

                x = line['x'][-1]  # last point
                line['x_dist'] = x2 - x
                points.append((line['x_dist'], (x+x2)/2, line['y']))
            else:
                line['x_dist'] = x2 - x1  # no points between
                points.append((line['x_dist'], (x1+x2)/2, line['y']))

        points.sort()
        return (points[-1][1], points[-1][2])

    def get_simple_point_inside(self, geom):
        mid_angle = middle_angle(self.min_angle, self.max_angle)
        c = (0.0, 0.0)
        p = point(c, self.max_dist + 5.0, mid_angle)
        line = Line(Element(start=c, end=p))
        points = []
        for e in self.area:
            points += e.intersect_line(line, rtol=geom.rtol, atol=geom.atol, include_end=True)

        if len(points) < 2:
            logger.warning("WARNING: get_simple_point_inside() failed (%s in %s)",
                           len(points), self.identifier())
            p = point(c, self.min_dist + (self.max_dist - self.min_dist) / 2, mid_angle)
            logger.warning("WARNING: try simple point %s", p)
            return p

        assert(len(points) > 1)
        points.sort()
        d1 = distance(c, points[0])
        d2 = distance(c, points[1])
        p = point(c, (d1 + d2) / 2, mid_angle)
        return p

    def get_point_inside(self, geom):
        """return point inside area"""
        y = (self.min_y + self.max_y) / 2
        p1 = (self.min_x - 5, y)
        p2 = (self.max_x + 5, y)
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
            if self.is_air():
                return self.get_best_point_inside(geom)
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

        if self.is_air():
            return self.get_best_point_inside(geom)
        return p_inside

    def render(self, renderer, color='black', with_nodes=False, fill=True):
        if fill:
            if self.render_fill(renderer):
                color = 'black'

        for e in self.area:
            e.render(renderer, color, with_nodes)
        return

    def render_fill(self, renderer):
        color = self.color()
        if not color:
            return False
        alpha = self.color_alpha()

        if self.is_circle():
            e = self.area[0]
            renderer.fill_circle(e.center, e.radius, color, alpha)
        else:
            nodes = [n for n in self.virtual_nodes(render=True)]
            x = [n[0] for n in nodes]
            y = [n[1] for n in nodes]
            renderer.fill(x, y, color, alpha)
        return True

    def magnet_arrow_length(self):
        if self.is_type(TYPE_MAGNET_AIRGAP):
            return (self.max_dist - self.min_dist) * 0.9
        if self.is_type(TYPE_MAGNET_RECT):
            return self.mag_width
        return 0.0

    def render_magnet_phi(self, renderer, length):
        if not self.is_magnet():
            return
        p1 = None
        p2 = None
        if self.is_type(TYPE_MAGNET_AIRGAP):
            mid = middle_angle(self.min_angle, self.max_angle)
            d = self.min_dist + (self.max_dist - self.min_dist) * 0.3
            p1 = point((0.0, 0.0), d, mid)
        if self.is_type(TYPE_MAGNET_RECT):
            x = self.min_x + (self.max_x - self.min_x) / 2
            y = self.min_y + (self.max_y - self.min_y) / 2
            p1 = (x, y)
        if p1 is None:
            return
        p2 = point(p1, length, self.phi)
        renderer.arrow(p1, p2, linewidth=1.2)
        return

    def render_legend(self, renderer):
        return renderer.new_legend_handle(self.color(),
                                          self.color_alpha(),
                                          self.legend())

    def remove_edges(self, g, ndec):
        r = 0
        for e in self.area:
            try:
                g.remove_edge(e.node1(ndec), e.node2(ndec))
                r+=1
            except Exception:
                continue
        return r

    def is_circle(self):
        e = self.area[0]
        if len(self.area) == 1:
            return is_Circle(e)

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

    def is_shaft_area(self, center):
        logger.debug("Begin of check shaft")

        #if not self.is_touching_both_sides():
        #    logger.debug("End of check shaft: don't touch both sides")
        #    return False

        if np.isclose(0.0, self.min_dist, rtol=1e-6, atol=1e-4):
            logger.debug("End of check shaft: ok (node in center)")
            return True

        for n in self.list_of_nodes():
            a = alpha_line(center, n)
            if np.isclose(self.min_angle, a):
                continue
            if np.isclose(self.max_angle, a):
                continue
            d = distance(center, n)
            if np.isclose(d, self.min_dist, atol=0.05):
                continue
            if np.isclose(d, self.max_dist, atol=0.05):
                continue
            logger.debug("End of check shaft: no")
            return False

        logger.debug("End of check shaft: ok")
        return True

    def get_magnet_line_angles(self):
        lines = [e for e in self.area if is_Line(e)]
        if len(lines) < 4:
            logger.debug("get_magnet_line_angles: only %s lines", len(lines))
            return []

        angles = []
        prev_angle = lines[0].get_positive_angle()
        logger.debug("first angle = %s", prev_angle)
        prev_length = lines[0].length()

        for line in lines[1:]:
            this_angle = line.get_positive_angle()
            logger.debug("next angle = %s", this_angle)
            this_length = line.length()

            if np.isclose(prev_angle, this_angle, rtol=1e-04, atol=1e-02):
                # same direction
                prev_length += this_length
            else:
                angles.append((prev_length, prev_angle))
                prev_angle = this_angle
                prev_length = this_length

        if not angles:
            logger.debug("get_magnet_line_angles: only one angle")
            return []

        this_length, this_angle = angles[0]
        if not np.isclose(prev_angle, this_angle, rtol=1e-04, atol=1e-02):
            angles.append((prev_length, prev_angle))
        else:
            prev_length += this_length
            angles[0] = (prev_length, prev_angle)

        l, first_angle = angles[0]
        l, last_angle = angles[-1]
        if np.isclose(first_angle, last_angle, rtol=1e-04, atol=1e-02):
            del angles[-1]
        return angles

    def get_magnet_phi(self, angles):
        if not angles:
            return 0.0

        angles.sort(reverse=True)
        # calculate orientation (no rectangle check)

        l, alpha = angles[0]
        phi = normalise_angle(alpha + np.pi/2)
        logger.debug("alpha = %s, phi = %s", alpha, phi)

        mid = middle_angle(self.min_angle, self.max_angle)
        angle = alpha_angle(mid, phi)
        logger.debug("phi=%s, mid=%s, angle=%s", phi, mid, angle)

        if np.isclose(mid, alpha, rtol=1e-3, atol=1e-3):
            phi = mid
            logger.debug("correction of phi=%s", phi)
        elif greater(angle, np.pi * 0.5, rtol=1e-5) and \
           less(angle, np.pi * 1.5, rtol=1e-5):
            phi = normalise_angle(phi + np.pi)
            logger.debug("correction of phi=%s", phi)

        logger.debug("phi of magnet %s is %s", self.identifier(), phi)
        return phi

    def get_magnet_orientation(self):
        logger.debug("get magnet orientation for %s", self.identifier())
        if self.is_type(TYPE_MAGNET_RECT):
            angles = self.get_magnet_line_angles()
            return self.get_magnet_phi(angles)

        if self.is_type(TYPE_MAGNET_AIRGAP):
            if self.close_to_endangle:
                if self.close_to_startangle:
                    return middle_angle(self.min_angle, self.max_angle)
                else:
                    return self.max_angle
            return middle_angle(self.min_angle, self.max_angle)

        return 0.0

    def is_magnet_rectangle(self):
        angles = self.get_magnet_line_angles()

        if len(angles) != 4:
            logger.debug("is_magnet_rectangle: %s angles, not 4", len(angles))
            return False

        for l, a in angles:
            logger.debug("+ magnet_rectangle: alpha=%s,  length=%s", a, l)

        length_0, angle_0 = angles[0]
        length_1, angle_1 = angles[1]
        length_2, angle_2 = angles[2]
        length_3, angle_3 = angles[3]

        if not np.isclose(angle_0, angle_2, rtol=1e-03, atol=0.05):
            logger.debug("is_magnet_rectangle: angles %s and %s not equal",
                         angle_0, angle_2)
            return False

        if not np.isclose(angle_1, angle_3, rtol=1e-03, atol=0.05):
            logger.debug("is_magnet_rectangle: angles %s and %s not equal",
                         angle_1, angle_3)
            return False

        if angle_0 > angle_1:
            a0 = angle_0
            a1 = angle_1 + np.pi/2
        else:
            a0 = angle_1
            a1 = angle_0 + np.pi/2
        if not np.isclose(a0, a1, rtol=1e-03, atol=0.05):
            logger.debug("is_magnet_rectangle: not a rectange (%s neq %s)", a0, a1)
            return False

        turn_left = False
        turn_right = False
        for n1, n2, e in self.list_of_elements():
            if is_Arc(e):
                if e.get_node_number(n1) == 1:
                    turn_left = True
                else:
                    turn_right = True

        if turn_left and turn_right:
            logger.debug("is_magnet_rectangle: arcs with different directions")
            return False

        self.phi = self.get_magnet_phi(angles)
        angles.sort()
        l, alpha = angles[0]
        self.mag_width = l
        logger.debug("Area %s is a rectangle with phi %s",
                     self.identifier(), self.phi)
        return True

    def around_windings(self, areas, geom):
        for a in areas:
            if a.is_winding():
                if not self.is_identical(a):
                    if self.is_inside(a, geom):
                        return True
                    elif self.is_touching(a):
                        return True
        return False

    def is_touching_areas(self, areas):
        for a in areas:
            if self.is_touching(a):
                return True
        return False

    def has_iron_separator(self):
        for e in self.area:
            if e.has_attribute('iron_sep'):
                return True
        return False

    def is_close_to_border(self, angle, border_angle):
        return np.isclose(angle, border_angle,
                          rtol=1e-03, atol=1e-03)

    def set_subregion_parameters(self,
                                 is_inner,
                                 min_radius,
                                 max_radius,
                                 startangle,
                                 endangle):
        ag_delta = (max_radius - min_radius) / 500.0
        if is_inner:
            self.close_to_ag = greater_equal(self.max_dist + ag_delta, max_radius)
            close_to_opposition = np.isclose(min_radius, self.min_dist)
            airgap_radius = max_radius
            opposite_radius = min_radius
        else:
            self.close_to_ag = less_equal(self.min_dist - ag_delta, min_radius)
            close_to_opposition = np.isclose(max_radius, self.max_dist)
            airgap_radius = min_radius
            opposite_radius = max_radius

        airgap_toleranz = (self.max_dist - self.min_dist) / 50.0  # 2%

        self.close_to_startangle = np.isclose(self.min_angle, startangle,
                                              1e-04, 1e-04)
        self.close_to_endangle = np.isclose(self.max_angle, endangle,
                                            1e-04, 1e-04)
        self.surface = self.area_size()

    def mark_stator_subregions(self,
                               is_inner,
                               stator_size,
                               mirrored,
                               alpha,
                               center,
                               r_in,
                               r_out):
        alpha = round(alpha, 6)

        if self.is_circle():
            self.type = TYPE_AIR  # air
            return self.type

        ag_delta = (r_out - r_in) / 500.0
        if is_inner:
            self.close_to_ag = greater_equal(self.max_dist + ag_delta, r_out)
            close_to_opposition = np.isclose(r_in, self.min_dist)
            airgap_radius = r_out
            opposite_radius = r_in
            airgap_toleranz = -(self.max_dist - self.min_dist) / 50.0  # 2%
        else:
            self.close_to_ag = less_equal(self.min_dist - ag_delta, r_in)
            close_to_opposition = np.isclose(r_out, self.max_dist)
            airgap_radius = r_in
            opposite_radius = r_out
            airgap_toleranz = (self.max_dist - self.min_dist) / 50.0  # 2%

        self.close_to_startangle = np.isclose(self.min_angle, 0.0,
                                              1e-04, 1e-04)
        self.close_to_endangle = np.isclose(self.max_angle, alpha,
                                            1e-04, 1e-04)
        self.surface = self.area_size()

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
        logger.debug(" - surface size       : %3.12f", self.surface)

        if self.has_iron_separator():
            logger.debug("***** iron (has iron separator)\n")
            self.type = TYPE_IRON  # iron
            return self.type

        if is_inner:
            # looking for shaft
            if close_to_opposition and not self.close_to_ag:
                if self.is_shaft_area(center):
                    self.type = TYPE_SHAFT  # shaft
                    logger.debug("***** shaft (close to opposition)\n")
                    return self.type

        if close_to_opposition:
            self.type = TYPE_YOKE  # iron yoke (Joch)
            logger.debug("***** iron yoke #1\n")
            return self.type

        if self.close_to_startangle and self.close_to_endangle:
            self.type = TYPE_YOKE  # iron yoke (Joch)
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
                self.type = TYPE_WINDINGS  # windings
            else:
                self.type = TYPE_AIR_OR_IRON  # air or iron near windings and near airgap?
                logger.debug("***** air or iron ??\n")
            return self.type

        if self.close_to_startangle:
            if self.is_half_circle(center, self.min_angle):
                self.type = TYPE_AIR  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        if self.close_to_endangle:
            if self.is_half_circle(center, self.max_angle):
                self.type = TYPE_AIR  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        def bad_winding_position():
            if is_inner:
                radius_third = airgap_radius - (airgap_radius - opposite_radius) * 0.33
                if self.max_dist < radius_third:
                    return True
            else:  # outer
                radius_third = airgap_radius + (opposite_radius - airgap_radius) * 0.33
                if self.min_dist > radius_third:
                    return True
            return False

        if self.min_angle > 0.001:
            if self.max_angle < alpha - 0.001:
                if bad_winding_position():
                    self.type = TYPE_WINDINGS_OR_AIR  # windings or air
                    logger.debug("***** windings or air #1\n")
                else:
                    self.type = TYPE_WINDINGS  # windings
                    logger.debug("***** windings #1\n")
                return self.type
            if mirrored:
                if bad_winding_position():
                    self.type = TYPE_WINDINGS_OR_AIR  # windings or air
                    logger.debug("***** windings or air #2\n")
                else:
                    self.type = TYPE_WINDINGS  # windings
                    logger.debug("***** windings #2\n")
                return self.type

            self.type = TYPE_AIR  # air
            logger.debug("***** air #3")

        if self.close_to_startangle or self.close_to_endangle:
            f = self.surface / stator_size
            if f < 0.02:  # area_size less then 2 percent of stator size
                # Luftloch
                self.type = TYPE_AIR  # air
                logger.debug("***** small area => air\n")
            else:
                self.type = TYPE_AIR_OR_IRON  # air or iron near windings and near airgap?
                logger.debug("***** air or iron close to border\n")
            return self.type

        logger.debug("***** air #4\n")
        return 0

    def mark_EESM_rotor_subregions(self,
                                   is_inner,
                                   mirrored,
                                   alpha,
                                   center,
                                   r_in,
                                   r_out,
                                   startangle,
                                   endangle):
        logger.debug("mark_EESM_rotor_subregions")

        alpha = round(alpha, 6)

        if self.is_circle():
            self.type = TYPE_AIR  # air
            logger.debug(">>> air is a circle")
            return self.type

        if is_inner:
            self.close_to_ag = np.isclose(r_out, self.max_dist, rtol=1e-9, atol=0.005)
            close_to_opposition = greater_equal(r_in * 1.05, self.min_dist)
            airgap_radius = r_out
            opposite_radius = r_in
            airgap_toleranz = -(self.max_dist - self.min_dist) / 50.0  # 2%
        else:
            self.close_to_ag = np.isclose(r_in, self.min_dist, rtol=1e-9, atol=0.005)
            close_to_opposition = greater_equal(self.max_dist * 1.05, r_out)
            airgap_radius = r_in
            opposite_radius = r_out
            airgap_toleranz = (self.max_dist - self.min_dist) / 50.0  # 2%

        self.close_to_startangle = np.isclose(self.min_angle, startangle,
                                              1e-04, 1e-04)
        self.close_to_endangle = np.isclose(self.max_angle, endangle,
                                            1e-04, 1e-04)

        logger.debug("\n***** mark_EESM_rotor_subregions [{}] *****"
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

        if self.has_iron_separator():
            logger.debug("***** iron (has iron separator)\n")
            self.type = TYPE_IRON  # iron
            return self.type

        if is_inner:
            # looking for shaft
            if close_to_opposition and not self.close_to_ag:
                logger.debug("-- check for shaft")
                if self.is_shaft_area(center):
                    self.type = TYPE_SHAFT  # shaft
                    logger.debug("***** shaft (close to opposition)\n")
                    return self.type

        if close_to_opposition or self.close_to_ag:
            self.type = TYPE_IRON  # iron
            logger.debug("***** iron (close to opposition)\n")
            return self.type

        self.type = TYPE_WINDINGS_OR_IRON  # windings or iron
        logger.debug("***** air (somewhere)\n")
        return self.type

    def mark_PMSM_rotor_subregions(self,
                                   is_inner,
                                   mirrored,
                                   alpha,
                                   center,
                                   r_in,
                                   r_out,
                                   startangle,
                                   endangle):
        logger.debug("mark_PMSM_rotor_subregions")

        alpha = round(alpha, 6)

        if self.is_circle():
            self.type = TYPE_AIR  # air
            logger.debug(">>> air is a circle")
            return self.type

        if is_inner:
            self.close_to_ag = np.isclose(r_out, self.max_dist, rtol=1e-9, atol=0.005)
            close_to_opposition = greater_equal(r_in * 1.05, self.min_dist)
            airgap_radius = r_out
            opposite_radius = r_in
            airgap_toleranz = -(self.max_dist - self.min_dist) / 50.0  # 2%
        else:
            self.close_to_ag = np.isclose(r_in, self.min_dist, rtol=1e-9, atol=0.005)
            close_to_opposition = greater_equal(self.max_dist * 1.05, r_out)
            airgap_radius = r_in
            opposite_radius = r_out
            airgap_toleranz = (self.max_dist - self.min_dist) / 50.0  # 2%

        self.close_to_startangle = np.isclose(self.min_angle, startangle,
                                              1e-04, 1e-04)
        self.close_to_endangle = np.isclose(self.max_angle, endangle,
                                            1e-04, 1e-04)

        logger.debug("\n***** mark_PMSM_rotor_subregions [{}] *****"
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

        if self.has_iron_separator():
            logger.debug("***** iron (has iron separator)\n")
            self.type = TYPE_IRON  # iron
            return self.type

        if is_inner:
            # looking for shaft
            if close_to_opposition and not self.close_to_ag:
                logger.debug("-- check for shaft")
                if self.is_shaft_area(center):
                    self.type = TYPE_SHAFT  # shaft
                    logger.debug("***** shaft (close to opposition)\n")
                    return self.type

        if close_to_opposition:
            self.type = TYPE_IRON  # iron
            logger.debug("***** iron (close to opposition)\n")
            return self.type

        if self.close_to_startangle and self.close_to_endangle:
            self.type = TYPE_IRON  # iron
            logger.debug("***** iron (close to both sides)\n")
            return self.type

        self.mag_rectangle = self.is_magnet_rectangle()

        if self.close_to_ag:
            mm = self.minmax_angle_dist_from_center(center,
                                                    airgap_radius +
                                                    airgap_toleranz)
            air_alpha = round(alpha_angle(mm[0], mm[1]), 3)
            logger.debug(" - air_alpha          : {}".format(air_alpha))

            if self.mag_rectangle:
                self.type = TYPE_MAGNET_RECT_NEAR_AIRGAP  # magnet near airgap
                logger.debug("***** magnet (airgap, embedded, phi={})\n".
                             format(self.phi))
                return self.type

            if air_alpha / alpha < 0.2:
                self.type = TYPE_MAGNET_OR_AIR  # air or magnet ?
                logger.debug("***** air #1 (close to airgap)\n")
                return self.type

            if air_alpha / alpha > 0.6:
                self.type = TYPE_MAGNET_AIRGAP  # magnet (no rectangle)
                logger.debug("***** magnet (close to airgap)\n")
            else:
                self.type = TYPE_MAGNET_OR_IRON  # iron or magnet ?
                logger.debug("***** iron or magnet(close to airgap)\n")
            return self.type

        if self.mag_rectangle:
            # phi is already calculated and set
            self.type = TYPE_MAGNET_RECT  # magnet embedded
            logger.debug("***** magnet (embedded, phi={})\n".format(
                self.phi))
            return self.type

        if not (self.close_to_startangle or self.close_to_endangle):
            self.type = TYPE_AIR  # air
            logger.debug("***** air (somewhere)\n")
            return self.type

        if self.close_to_startangle:
            if self.is_half_circle(center, self.min_angle):
                self.type = TYPE_AIR  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        if self.close_to_endangle:
            if self.is_half_circle(center, self.max_angle):
                self.type = TYPE_AIR  # air
                logger.debug("***** air (part of a circle)\n")
                return self.type

        self.type = TYPE_AIR  # air
        logger.debug("***** air (remains)\n")
        return self.type

    def mark_unknown_subregions(self, mirrored, alpha,
                                center, r_in, r_out):
        logger.debug("mark_unknown_subregions")

        if self.is_circle():
            self.type = TYPE_AIR  # air
            logger.debug(">>> air is a circle")
            return self.type

        self.close_to_startangle = np.isclose(self.min_angle, 0.0)
        self.close_to_endangle = np.isclose(self.max_angle, alpha)

        if self.is_magnet_rectangle():
            self.type = TYPE_MAGNET_RECT  # magnet embedded
            logger.debug(">>> magnet embedded")
            return self.type

        close_to_max_radius = np.isclose(r_out, self.max_dist)
        close_to_min_radius = np.isclose(r_in, self.min_dist)

        if close_to_max_radius and close_to_min_radius:
            self.type = TYPE_IRON  # iron
            logger.debug(">>> iron close to min- and max-radius")
            return self.type

        if self.close_to_startangle and self.close_to_endangle:
            self.type = TYPE_IRON  # iron
            logger.debug(">>> iron close to start- and end-angle")
            return self.type

        self.type = TYPE_AIR  # air
        logger.debug(">>> air remains")
        return self.type

    def mark_airgap_corners(self, start_cp, end_cp):
        for n in self.list_of_nodes():
            if self.close_to_startangle:
                if points_are_close(n, start_cp):
                    self.close_to_ag_startcorner = True
            if self.close_to_endangle:
                if points_are_close(n, end_cp):
                    self.close_to_ag_endcorner = True

    def area_size(self):
        if len(self.area) == 0:
            return 0.0
        if self.number_of_elements() < 2:
            e = self.area[0]
            if not is_Circle(e):
                return 0.0
            return np.pi * e.radius**2

        nodes = [n for n in self.list_of_nodes()]
        return area_size(nodes)

    def set_surface(self, mirrored):
        self.surface = self.area_size()
        if self.close_to_endangle and mirrored:
            self.surface = self.area_size() * 2.0
        else:
            self.surface = self.area_size()

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

    def nested_areas_inside(self):
        for id, a in self.areas_inside.items():
            yield id
            for i in a.nested_areas_inside():
                yield i

    def list_of_nested_areas_inside(self):
        for id, a in self.areas_inside.items():
            for i in a.nested_areas_inside():
                yield i

    def mirror_area(self, center, axis_m, axis_n):
        elements = []
        for e in self.area:
            el = e.mirror_shape(center, axis_m, axis_n)
            if el:
                elements.append(el)
        area = Area(elements, center, 0.0)
        area.type = self.type
        return area

    def __str__(self):
        return "Area {}\n".format(self.id) + \
            "distance...............: from {} to {}\n".\
            format(round(self.min_dist, 4), round(self.max_dist, 4)) + \
            "height.................: {}\n".format(self.height) + \
            "alpha..................: {}\n".format(self.alpha) + \
            "angle..................: from {} to {}\n".\
            format(round(self.min_angle, 6), round(self.max_angle, 6)) + \
            "delta..................: {}\n".format(self.delta) + \
            "number.................: {}\n".format(self.count) + \
            "equal areas............: {}\n".format(len(self.equal_areas)) + \
            "symmetry...............: {}\n".format(self.symmetry) + \
            "symmetry type..........: {}\n".format(self.sym_type) + \
            "close to airgap........: {}\n".format(self.close_to_ag) + \
            "close to startangle....: {}\n".format(self.close_to_startangle) + \
            "close to endangle......: {}\n".format(self.close_to_endangle) + \
            "close to ag startcorner: {}\n".format(self.close_to_ag_startcorner) + \
            "close to ag endcorner..: {}\n".format(self.close_to_ag_endcorner)
