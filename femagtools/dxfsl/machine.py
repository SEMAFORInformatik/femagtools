# -*- coding: utf-8 -*-
"""
    femagtools.dxfsl.machine
    ~~~~~~~~~~~~~~~~~~~~~~~~

    a machine consists of 2 parts and has a geometry 

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function
import numpy as np
import logging
from .shape import Element, Circle, Arc, Line, Shape
from .corner import Corner
from .functions import point, points_are_close, distance
from .functions import alpha_angle, normalise_angle, middle_angle, third_angle
from .functions import line_m, line_n, mirror_point
from .functions import within_interval, part_of_circle
from .functions import less, less_equal, greater, greater_equal
logger = logging.getLogger('femagtools.geom')


#############################
#          Machine          #
#############################

class Machine(object):
    def __init__(self, geom, center, radius, startangle=0, endangle=0):
        self.geom = geom
        self.center = center
        self.radius = radius
        self.startangle = startangle
        self.endangle = endangle
        self.mirror_orig_geom = None
        self.mirror_geom = None
        self.mirror_startangle = 0.0
        self.mirror_endangle = 0.0
        self.part = self.part_of_circle()
        self.airgaps = []
        self.airgap_radius = 0.0
        self.airgap2_radius = 0.0
        self.geom.center = center
        self.previous_machine = None

    def __str__(self):
        return "Machine\n" + \
               "Center: ({})\n".format(self.center) + \
               "Radius: {}\n".format(self.radius) + \
               "Angles: Start={}, End={}\n".format(self.startangle,
                                                   self.endangle) + \
               "Mirror: {}\n".format(self.mirror_geom is not None)

    def is_a_machine(self):
        return self.radius > 0.0

    def is_in_middle(self):
        return self.radius > 0.0 and \
            points_are_close(self.center, [0.0, 0.0], 1e-8)

    def is_full(self):
        return self.radius > 0.0 and \
            self.startangle == 0.0 and self.endangle == 0.0

    def is_half_up(self):
        return self.radius > 0.0 and \
               np.isclose(self.startangle, 0.0, 1e-8) and \
               (np.isclose(self.endangle, np.pi, 1e-8) or
                np.isclose(self.endangle, -np.pi, 1e-8))

    def is_half_down(self):
        return self.radius > 0.0 and \
               (np.isclose(self.endangle, np.pi, 1e-8) or
                np.isclose(self.endangle, -np.pi, 1e-8)) and \
               np.isclose(self.startangle, 0.0, 1e-8)

    def is_half_left(self):
        return self.radius > 0.0 and \
               np.isclose(self.startangle, np.pi/2, 1e-8) and \
               np.isclose(self.endangle, -np.pi/2, 1e-8)

    def is_half_right(self):
        return self.radius > 0.0 and \
               np.isclose(self.startangle, -np.pi/2, 1e-8) and \
               np.isclose(self.endangle, np.pi/2, 1e-8)

    def is_half(self):
        return self.is_half_up() or self.is_half_down() or \
               self.is_half_left() or self.is_half_right()

    def is_quarter(self):
        return self.radius > 0.0 and \
               np.isclose(alpha_angle(self.startangle, self.endangle), np.pi/2)

    def is_startangle_zero(self):
        return np.isclose(self.startangle, 0.0)

    def move_to_middle(self):
        if not self.is_in_middle():
            if self.radius > 0.0:
                offset = [-(self.center[0]), -(self.center[1])]
                self.geom.move(offset)
                self.set_center(0.0, 0.0)
                self.geom.clear_cut_lines()
                return True
        return False

    def set_center(self, x, y):
        self.center[0] = x
        self.center[1] = y
        self.geom.center = [x, y]

    def set_radius(self, radius):
        self.radius = radius

    def set_kind(self, kind):
        self.geom.kind = kind

    def set_inner(self):
        self.geom.is_inner = True

    def set_outer(self):
        self.geom.is_outer = True

    def clear_cut_lines(self):
        self.geom.clear_cut_lines()
        if self.mirror_geom is not None:
            self.mirror_geom.clear_cut_lines()

    def cut_is_possible(self, r_in, r_out):
        if r_in > 0.0:
            if less(self.geom.min_radius, r_in, rtol=0.0001):
                if less_equal(self.geom.max_radius, r_in, rtol=0.0001):
                    return False
                return True
        if r_out > 0.0:
            if greater(self.geom.max_radius, r_out, rtol=0.0001):
                if greater_equal(self.geom.min_radius, r_out, rtol=0.0001):
                    return False
                return True
        return False

    def cut(self, r_in, r_out):
        radius_in = r_in
        radius_out = r_out
        if r_out == 0.0:
            radius_out = self.radius + 10

        clone = self.geom.copy_shape(self.center,
                                     self.radius,
                                     0.0,
                                     2*np.pi,
                                     radius_in,
                                     radius_out,
                                     False,
                                     append_inner=(r_in > 0.0),
                                     append_outer=(r_out > 0.0))

        if r_out == 0.0:
            r_out = self.radius

        m = Machine(clone, self.center, r_out,
                    self.startangle, self.endangle)

        m.mirror_geom = self.mirror_geom
        m.part = self.part

        m.set_minmax_radius()
        m.create_auxiliary_lines()
        m.set_alfa_and_corners()
        m.set_kind(self.geom.kind)
        return m

    def copy(self, startangle, endangle,
             airgap=False, inside=True, split=False):
        if airgap and self.airgap_radius > 0.0:
            if inside:
                if self.airgap2_radius > 0.0:
                    new_radius = min(self.airgap_radius, self.airgap2_radius)
                else:
                    new_radius = self.airgap_radius
                clone = self.geom.copy_shape(self.center,
                                             self.radius,
                                             startangle, endangle,
                                             0.0, new_radius, split)
            else:
                new_radius = self.radius
                gap_radius = max(self.airgap_radius, self.airgap2_radius)
                clone = self.geom.copy_shape(self.center,
                                             self.radius,
                                             startangle, endangle,
                                             gap_radius, self.radius+9999,
                                             split)

            circ = Circle(Element(center=self.center,
                                  radius=self.airgap_radius))
            clone.add_cut_line(circ)
        else:
            new_radius = self.radius
            clone = self.geom.copy_shape(self.center, self.radius,
                                         startangle, endangle, 0.0,
                                         self.radius+9999, split)

        if not np.isclose(normalise_angle(startangle),
                          normalise_angle(endangle), 0.0):
            start_p = point(self.center, self.radius+5, startangle)
            start_line = Line(Element(start=self.center, end=start_p))
            clone.add_cut_line(start_line)

            end_p = point(self.center, self.radius+5, endangle)
            end_line = Line(Element(start=self.center, end=end_p))
            clone.add_cut_line(end_line)

        if not np.isclose(alpha_angle(startangle, endangle), 2*np.pi):
            return Machine(clone, self.center, new_radius,
                           startangle, endangle)
        else:
            # Der Originalwinkel bleibt bestehen
            return Machine(clone, self.center, new_radius,
                           self.startangle, self.endangle)

    def full_copy(self):
        clone = self.geom.copy_shape(self.center, self.radius,
                                     0.0, 2*np.pi,
                                     0.0, self.radius+9999)
        return clone.get_machine()

    def copy_mirror(self, startangle, midangle, endangle):
        geom1 = self.geom.copy_shape(self.center,
                                     self.radius,
                                     startangle,
                                     midangle,
                                     0.0,
                                     self.radius+9999,
                                     rtol=1e-08,
                                     atol=1e-08)
        geom2 = self.geom.copy_shape(self.center,
                                     self.radius,
                                     midangle,
                                     endangle,
                                     0.0,
                                     self.radius+9999,
                                     rtol=1e-08,
                                     atol=1e-08)

        machine = Machine(geom1, self.center, self.radius,
                          startangle, midangle)
        machine.mirror_orig_geom = self.geom
        machine.mirror_geom = geom2
        machine.mirror_geom.center = self.center
        machine.mirror_startangle = midangle
        machine.mirror_endangle = endangle
        return machine

    def has_mirrored_windings(self):
        if not self.is_mirrored():
            return False
        return self.geom.area_close_to_endangle(2) > 0

    def undo_mirror(self):
        assert(self.is_mirrored())
        assert(self.previous_machine)
        self.previous_machine.set_minmax_radius()
        # self.previous_machine.complete_hull()
        self.set_alfa_and_corners()
        self.previous_machine.create_auxiliary_lines()
        self.previous_machine.set_kind(self.geom.kind)
        return self.previous_machine

    def rotate_to(self, new_startangle):
        if np.isclose(new_startangle, self.startangle):
            return

        if points_are_close(self.center, [0.0, 0.0]):
            angle = new_startangle - self.startangle
            self.geom.rotate(angle)
            self.startangle = new_startangle
            self.endangle += angle

    def airgap(self, correct_airgap=0.0, correct_airgap2=0.0, atol=0.1):
        logger.debug('locking for airgap')
        self.airgap_radius = 0.0
        self.airgap2_radius = 0.0

        if np.isclose(self.radius, 0.0):
            logger.debug('no radius')
            return False

        if correct_airgap < 0:
            logger.debug('no airgap')
            return False  # no airgap

        self.airgaps = []
        airgaps = self.geom.detect_airgaps(self.center,
                                           self.startangle,
                                           self.endangle, atol)

        alpha = alpha_angle(self.startangle, self.endangle)

        if len(airgaps) == 1:
            self.airgaps = airgaps
        elif len(airgaps) > 0:
            lower_radius = -1.0
            upper_radius = -1.0

            for g in airgaps:
                if np.isclose(g[0], upper_radius):
                    if not self.geom.delete_airgap_circle(self.center,
                                                          lower_radius,
                                                          upper_radius,
                                                          g[1],
                                                          alpha):
                        lower_radius = g[0]
                else:
                    if lower_radius > 0.0:
                        self.airgaps.append((lower_radius, upper_radius))
                    lower_radius = g[0]
                upper_radius = g[1]
            self.airgaps.append((lower_radius, upper_radius))

        if len(self.airgaps) > 0:
            airgap_candidates = []
            for g in self.airgaps:
                gap_radius = round((g[0]+g[1])/2.0, 6)
                gap_dist = g[1] - g[0]
                circle = Circle(Element(center=self.center,
                                        radius=gap_radius))
                ok, borders = self.geom.is_airgap(self.center,
                                                  self.radius,
                                                  self.startangle,
                                                  self.endangle,
                                                  circle, atol)
                if not ok:
                    logger.error("FATAL: No Airgap with radius {}".
                                 format(gap_radius))
                    print("FATAL: No Airgap with radius {}".
                          format(gap_radius))
                    self.geom.airgaps.append(circle)
                    return True  # bad exit

                airgap_candidates.append((borders, circle, gap_dist))
                self.geom.airgaps.append(circle)

                if correct_airgap > 0.0:
                    if within_interval(correct_airgap, g[0], g[1], 0.0, 0.0):
                        self.airgap_radius = gap_radius  # ok

                if correct_airgap2 > 0.0:
                    if within_interval(correct_airgap2, g[0], g[1], 0.0, 0.0):
                        self.airgap2_radius = gap_radius  # ok

        if correct_airgap > 0.0 and self.airgap_radius == 0.0:
            logger.error("No airgap with radius {} found"
                         .format(correct_airgap))
            self.show_airgap_candidates(airgap_candidates, False)
            return True  # bad exit

        if correct_airgap2 > 0.0 and self.airgap2_radius == 0.0:
            logger.error("No airgap2 with radius {} found"
                         .format(correct_airgap2))
            self.show_airgap_candidates(airgap_candidates, False)
            return True  # bad exit

        if len(self.airgaps) == 0:
            logger.debug('No airgap found')
            return False  # no airgaps found

        if self.airgap_radius > 0.0:
            return False  # correct airgap set

        gaps = [c for b, c, d in airgap_candidates if b == 0]

        if len(gaps) == 1:  # one candidate without border intersection
            self.airgap_radius = gaps[0].radius
            return False  # ok

        if len(airgap_candidates) == 1:  # one candidate found
            self.airgap_radius = airgap_candidates[0][1].radius
            return False  # ok

        self.airgap_radius = self.show_airgap_candidates(airgap_candidates,
                                                         True)
        return False  # ok

    def show_airgap_candidates(self, airgap_candidates, get_one):
        if get_one:
            logger.info("{} airgap candidate(s) found:"
                        .format(len(airgap_candidates)))
        else:
            print("{} airgap candidate(s) found:"
                  .format(len(airgap_candidates)))
        dist = 999
        circle = None
        pos_list = []
        for b, c, d in airgap_candidates:
            if get_one:
                logger.info(" --- {}   (width={})".format(c.radius, d))
            else:
                print(" --- {}   (width={})".format(c.radius, d))
            if d < dist:
                dist = d
                circle = c

            inner_pc = (c.radius - self.geom.min_radius) / \
                       (self.geom.max_radius - self.geom.min_radius)
            pos = np.abs(inner_pc * 100 - 50)
            logger.debug("Abstand Mitte = {} %".format(pos))
            if pos < 20:
                pos_list.append([pos, d, c])

        if get_one:
            if pos_list:
                dist_list = [[d, c] for pos, d, c in pos_list]
                dist_list.sort()
                circle = dist_list[0][1]

            logger.info("airgap {} prefered".format(circle.radius))
            return circle.radius

        print("Use options --airgap/--airgap2 <float> to specify")
        return None

    def has_airgap(self):
        return self.airgap_radius > 0.0

    def airgap_x(self):
        return self.airgap_radius

    def airgap_y(self):
        return 0.1

    def part_of_circle(self, pos=3):
        return part_of_circle(self.startangle, self.endangle, pos)

    def delete_center_circle(self):
        gaps = self.geom.get_gaplist(self.center)
        if len(gaps) < 2:
            return

        first_gap = gaps[0][1]
        second_gap = gaps[1][0]
        if first_gap != second_gap:
            return

        first_dist = gaps[0][1]
        second_dist = gaps[1][1] - gaps[1][0]
        if first_dist < 1.0:
            if second_dist / first_dist > 10.0:
                self.geom.delete_circle((0.0, 0.0), first_dist)

    def repair_hull(self):
        logger.debug('repair_hull')
        if self.is_full() and not self.has_airgap():
            self.delete_center_circle()

        if self.startangle == self.endangle:
            logger.info('end of repair_hull: circle')
            return

        self.repair_hull_geom(self.geom, self.startangle, self.endangle)

        if self.mirror_geom:
            self.repair_hull_geom(self.mirror_geom,
                                  self.mirror_startangle,
                                  self.mirror_endangle)
        logger.debug('end of repair_hull')

    def repair_hull_geom(self, geom, startangle, endangle):
        logger.debug('repair_hull_geom')

        c_corner = Corner(self.center, self.center)
        start_corners = geom.get_corner_list(self.center, startangle)
        end_corners = geom.get_corner_list(self.center, endangle)

        geom.repair_hull_line(self.center,
                              startangle, start_corners,
                              c_corner in end_corners)
        geom.repair_hull_line(self.center,
                              endangle, end_corners,
                              c_corner in start_corners)
        logger.debug('end of repair_hull_geom')

    def set_minmax_radius(self):
        self.geom.set_minmax_radius(self.center)

    def complete_hull(self, is_inner, is_outer):
        logger.info('complete_hull')
        start_corners = self.geom.complete_hull_line(self.center,
                                                     self.startangle)
        end_corners = self.geom.complete_hull_line(self.center,
                                                   self.endangle)

        if start_corners[0].is_new_point or end_corners[0].is_new_point:
            self.geom.complete_hull_arc(self.center,
                                        self.startangle, start_corners[0],
                                        self.endangle, end_corners[0],
                                        self.geom.min_radius)

        if start_corners[1].is_new_point or end_corners[1].is_new_point:
            self.geom.complete_hull_arc(self.center,
                                        self.startangle, start_corners[1],
                                        self.endangle, end_corners[1],
                                        self.geom.max_radius)

        self.set_alfa_and_corners()

    def create_auxiliary_lines(self):
        self.geom.create_auxiliary_lines(self.startangle, self.endangle)

    def set_alfa_and_corners(self):
        self.geom.start_corners = self.geom.get_corner_nodes(self.center,
                                                             self.startangle)
        self.geom.end_corners = self.geom.get_corner_nodes(self.center,
                                                           self.endangle)
        self.geom.alfa = alpha_angle(self.startangle, self.endangle)

        if self.mirror_geom is not None:
            self.geom.mirror_corners = self.geom.end_corners

    def is_mirrored(self):
        return self.mirror_geom is not None

    def num_of_layers(self):
        w = self.geom.num_of_windings()
        if w > 0 and self.is_mirrored():
            return w*2
        return w

    def find_symmetry(self, sym_tolerance):
        if self.radius <= 0.0:
            return False

        return self.geom.find_symmetry(self.center, self.radius,
                                       self.startangle, self.endangle,
                                       sym_tolerance)

    def get_symmetry_slice(self):
        logger.debug("begin get_symmetry_slice")
        if not self.geom.has_symmetry_area():
            logger.debug("end get_symmetry_slice: no symmetry area")
            return None

        machine_slice = self.copy(self.geom.symmetry_startangle(),
                                  self.geom.symmetry_endangle())
        machine_slice.clear_cut_lines()
        machine_slice.repair_hull()
        machine_slice.rotate_to(0.0)
        machine_slice.set_alfa_and_corners()
        logger.debug("end get_symmetry_slice: angle start: {}, end: {}"
                     .format(self.geom.symmetry_startangle(),
                             self.geom.symmetry_endangle()))
        return machine_slice

    def get_third_symmetry_mirror(self):
        logger.debug("begin get_third_symmetry_mirror")
        first_thirdangle = third_angle(self.startangle, self.endangle)
        second_thirdangle = middle_angle(first_thirdangle, self.endangle)

        machine_mirror_1 = self.copy_mirror(self.startangle,
                                            first_thirdangle,
                                            second_thirdangle)
        machine_mirror_1.clear_cut_lines()
        machine_mirror_1.repair_hull()
        machine_mirror_1.set_alfa_and_corners()
        if not machine_mirror_1.check_symmetry_graph(0.001, 0.05):
            logger.debug("end get_third_symmetry_mirror: no mirror first third")
            return None

        machine_mirror_2 = self.copy_mirror(first_thirdangle,
                                            second_thirdangle,
                                            self.endangle)
        machine_mirror_2.clear_cut_lines()
        machine_mirror_2.repair_hull()
        machine_mirror_2.set_alfa_and_corners()
        if not machine_mirror_2.check_symmetry_graph(0.001, 0.05):
            logger.debug("end get_third_symmetry_mirror: no mirror second third")
            return None

        machine_mirror_1.previous_machine = self
        logger.debug("end get_third_symmetry_mirror: ok")
        return machine_mirror_1

    def get_symmetry_mirror(self):
        logger.debug("begin get_symmetry_mirror")
        if self.part == 1:
            # a complete machine
            startangle = 0.0
            endangle = 0.0
            midangle = np.pi
        else:
            startangle = self.startangle
            endangle = self.endangle
            midangle = middle_angle(self.startangle, self.endangle)
            machine_mirror = self.get_third_symmetry_mirror()
            if machine_mirror:
                logger.debug("end get_symmetry_mirror: third found")
                return machine_mirror

        logger.debug(" - angles: start: {}, mid: {}, end: {}"
                    .format(startangle, midangle, endangle))
        machine_mirror = self.copy_mirror(startangle, midangle, endangle)
        machine_mirror.clear_cut_lines()
        machine_mirror.repair_hull()
        machine_mirror.set_alfa_and_corners()
        if machine_mirror.check_symmetry_graph(0.001, 0.05):
            machine_mirror.previous_machine = self
            machine_mirror.rotate_to(0.0)
            machine_mirror.set_alfa_and_corners()
            logger.debug("end get_symmetry_mirror: found")
            return machine_mirror

        logger.debug("end get_symmetry_mirror: no mirror")
        return None

    def get_symmetry_part(self):
        if self.is_mirrored():
            return self.part/2
        else:
            return self.part

    def check_symmetry_graph(self, rtol, atol):
        logger.debug("begin check_symmetry_graph")
        axis_p = point(self.center, self.radius, self.mirror_startangle)
        axis_m = line_m(self.center, axis_p)
        axis_n = line_n(self.center, axis_m)

        def is_node_available(n, nodes):
            mirror_p = mirror_point(n, self.center, axis_m, axis_n)
            for p in nodes:
                if points_are_close(p, mirror_p, rtol, atol):
                    return True
            return False

        def get_hit_factor(nodes1, nodes2):
            hit = 0
            if not nodes1:
                return 0.0

            for n in nodes1:
                if is_node_available(n, nodes2):
                    hit += 1
                else:
                    d = distance(self.center, n)
                    logger.debug(" -- r={}, d={}".format(self.radius, d))
                    if np.isclose(d, self.radius, rtol=0.075, atol=atol):
                        # very tolerant
                        logger.debug("NO HIT FOR {} ON OUTER HULL".format(n))
                        hit += 1

            return float(hit) / len(nodes1)

        hit_factor1 = get_hit_factor(self.geom.g.nodes(),
                                     self.mirror_geom.g.nodes())
        logger.debug("=> hit_factor1 = {}".format(hit_factor1))
        if hit_factor1 < 0.9:
            return False  # not ok

        hit_factor2 = get_hit_factor(self.mirror_geom.g.nodes(),
                                     self.geom.g.nodes())
        logger.debug("=> hit_factor2 = {}".format(hit_factor2))
        if hit_factor2 < 0.9:
            return False  # not ok

        if hit_factor1 < 0.93 and hit_factor2 < 0.93:
            return False  # not ok

        logger.debug("end check_symmetry_graph: ok")
        return True

    def sync_with_counterpart(self, cp_machine):
        self.geom.sym_counterpart = cp_machine.get_symmetry_part()
        self.geom.sym_part = self.get_symmetry_part()
        cp_machine.geom.sym_counterpart = self.get_symmetry_part()
        cp_machine.geom.sym_part = cp_machine.get_symmetry_part()

    def search_subregions(self):
        self.geom.search_subregions()

    def delete_tiny_elements(self, mindist):
        self.geom.delete_tiny_elements(mindist)

    def create_mirror_lines_outside_windings(self):
        if not self.geom.has_areas_touching_both_sides():
            return

        midangle = middle_angle(self.startangle, self.endangle)
        pts = self.geom.split_and_get_intersect_points(self.center,
                                                       self.radius+10,
                                                       midangle)
        pts.sort()
        if self.geom.create_lines_outside_windings(pts):
            self.geom.area_list = []
            logger.debug("create subregions again")
            self.geom.create_list_of_areas(crunch=True)
            self.geom.search_subregions()
