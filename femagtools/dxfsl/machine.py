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
import logging
import sys
from .shape import Element, Circle, Line
from .functions import point, points_are_close
from .functions import alpha_angle, normalise_angle, middle_angle
from .functions import line_m, line_n, mirror_point
from .functions import within_interval, part_of_circle
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
        geom1 = self.geom.copy_shape(self.center, self.radius,
                                     startangle, midangle,
                                     0.0, self.radius+9999)
        geom2 = self.geom.copy_shape(self.center, self.radius,
                                     midangle, endangle,
                                     0.0, self.radius+9999)
        machine = Machine(geom1, self.center, self.radius,
                          startangle, midangle)
        machine.mirror_orig_geom = self.geom
        machine.mirror_geom = geom2
        machine.mirror_geom.center = self.center
        machine.mirror_startangle = midangle
        machine.mirror_endangle = endangle
        return machine

    def undo_mirror(self):
        if self.is_mirrored():
            self.endangle = self.mirror_endangle
            self.mirror_orig_geom.min_radius = self.geom.min_radius
            self.mirror_orig_geom.max_radius = self.geom.max_radius
            self.mirror_orig_geom.kind = self.geom.kind
            self.geom = self.mirror_orig_geom
            self.mirror_orig_geom = None
            self.mirror_geom = None
            self.mirror_startangle = 0.0
            self.mirror_endangle = 0.0
            self.part = self.part_of_circle()
            self.set_alfa_and_corners()
            self.geom.create_list_of_areas()

    def rotate_to(self, new_startangle):
        if np.isclose(new_startangle, self.startangle):
            return

        if points_are_close(self.center, [0.0, 0.0]):
            angle = new_startangle - self.startangle
            self.geom.rotate(angle)
            self.startangle = new_startangle
            self.endangle += angle

    def airgap(self, correct_airgap=0.0, correct_airgap2=0.0, atol=0.05):
        self.airgap_radius = 0.0
        self.airgap2_radius = 0.0

        if np.isclose(self.radius, 0.0):
            return

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
            num_airgaps = 0
            for g in self.airgaps:
                gap_radius = round((g[0]+g[1])/2.0, 6)

                if correct_airgap == 0.0 or \
                   within_interval(correct_airgap, g[0], g[1], 0.0, 0.0):
                    circle = Circle(Element(center=self.center,
                                            radius=gap_radius))
                    if self.geom.is_airgap(self.center, self.radius,
                                           self.startangle,
                                           self.endangle, circle, atol):
                        self.geom.airgaps.append(circle)
                        num_airgaps += 1
                        self.airgap_radius = gap_radius
                    else:
                        logger.debug("DESASTER: No Airgap with radius {}".
                                     format(gap_radius))
                        print("DESASTER: No Airgap with radius {}".
                              format(gap_radius))
                        sys.exit(1)

                if correct_airgap2 > 0.0 and \
                   within_interval(correct_airgap2, g[0], g[1], 0.0, 0.0):
                    circle = Circle(Element(center=self.center,
                                            radius=gap_radius))
                    if self.geom.is_airgap(self.center, self.radius,
                                           self.startangle, self.endangle,
                                           circle, atol):
                        self.airgap2_radius = gap_radius
                    else:
                        logger.debug("DESASTER: No Airgap with radius {}".
                                     format(gap_radius))
                        print("DESASTER: No Airgap with radius {}".
                              format(gap_radius))
                        sys.exit(1)

            if num_airgaps == 1:
                return

            if num_airgaps > 1:
                airgaps = []
                for g in self.airgaps:
                    lower_circles = self.geom.get_circles(self.center, g[0])
                    upper_circles = self.geom.get_circles(self.center, g[1])
                    sum_lower_angle = self.geom.alpha_of_circles(lower_circles,
                                                                 self.center)
                    sum_upper_angle = self.geom.alpha_of_circles(upper_circles,
                                                                 self.center)
                    if(sum_lower_angle / alpha > 0.75 and
                       sum_upper_angle / alpha > 0.75):
                        airgaps.append(g)

                if len(airgaps) == 1:
                    g = airgaps[0]
                    self.airgap_radius = round((g[0]+g[1])/2.0, 6)
                    self.geom.airgaps = airgaps
                else:
                    print("More than one airgap candidate found:")
                    for c in self.geom.airgaps:
                        print(" --- {}".format(c.radius))
                        print("Use options --airgap/--airgap2 <float> to specify")
                    sys.exit(1)
            else:
                self.airgap_radius = 0.0

    def has_airgap(self):
        return self.airgap_radius > 0.0

    def airgap_x(self):
        return self.airgap_radius

    def airgap_y(self):
        return 0.1

    def part_of_circle(self, pos=3):
        return part_of_circle(self.startangle, self.endangle, pos)

    def repair_hull(self):
        self.geom.repair_hull_line(self.center, self.startangle)
        self.geom.repair_hull_line(self.center, self.endangle)

        if self.mirror_geom:
            self.mirror_geom.repair_hull_line(self.center,
                                              self.mirror_startangle)
            self.mirror_geom.repair_hull_line(self.center,
                                              self.mirror_endangle)

    def complete_hull(self):
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
        self.geom.create_auxiliary_lines(self.endangle)

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

    def find_symmetry(self, sym_tolerance):
        if self.radius <= 0.0:
            return False

        return self.geom.find_symmetry(self.center, self.radius,
                                       self.startangle, self.endangle,
                                       sym_tolerance)

    def get_symmetry_slice(self):
        if not self.geom.has_symmetry_area():
            return None

        machine_slice = self.copy(self.geom.symmetry_startangle(),
                                  self.geom.symmetry_endangle())
        machine_slice.clear_cut_lines()
        machine_slice.repair_hull()
        machine_slice.rotate_to(0.0)
        machine_slice.set_alfa_and_corners()
        return machine_slice

    def get_symmetry_mirror(self):
        if self.part == 1:
            # ein ganzer Motor
            startangle = 0.0
            endangle = 0.0
            midangle = np.pi
        else:
            startangle = self.startangle
            endangle = self.endangle
            midangle = middle_angle(self.startangle, self.endangle)

        machine_mirror = self.copy_mirror(startangle, midangle, endangle)
        machine_mirror.clear_cut_lines()
        machine_mirror.repair_hull()
        machine_mirror.set_alfa_and_corners()
        if machine_mirror.check_symmetry_graph(0.1, 0.1):
            return machine_mirror
        return None

    def get_symmetry_part(self):
        if self.is_mirrored():
            return self.part/2
        else:
            return self.part

    def check_symmetry_graph(self, rtol, atol):
        axis_p = point(self.center, self.radius, self.mirror_startangle)
        axis_m = line_m(self.center, axis_p)
        axis_n = line_n(self.center, axis_m)

        def is_node_available(n, nodes):
            mirror_p = mirror_point(n, self.center, axis_m, axis_n)
            for p in nodes:
                if points_are_close(p, mirror_p, rtol, atol):
                    return True
            return False

        hit = 0
        nodes = self.geom.g.nodes()
        for n in nodes:
            if is_node_available(n, self.mirror_geom.g.nodes()):
                hit += 1

        hit_factor = hit / len(nodes)
        ok = hit_factor > 0.9
#        print("Nodes = {}, Match={} => ok={}".format(len(nodes), hit, ok))
        return ok

    def sync_with_counterpart(self, cp_machine):
        self.geom.sym_counterpart = cp_machine.get_symmetry_part()
        self.geom.sym_part = self.get_symmetry_part()
        cp_machine.geom.sym_counterpart = self.get_symmetry_part()
        cp_machine.geom.sym_part = cp_machine.get_symmetry_part()

    def search_subregions(self):
        # print("search_subregions\n{}".format(self))
        self.geom.search_subregions()
        return
