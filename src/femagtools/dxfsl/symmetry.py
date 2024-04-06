# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.symmetry
  ~~~~~~~~~~~~~~~~~~~~~~~~~

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function

import numpy as np
import logging
import sys
from femagtools.dxfsl.shape import Element, Line
from femagtools.dxfsl.area import Area
from femagtools.dxfsl.functions import alpha_angle, positive_angle, is_same_angle
from femagtools.dxfsl.functions import min_angle, max_angle, gcd, point
from femagtools.dxfsl.functions import less_equal, less, points_are_close
from femagtools.dxfsl.functions import line_m, line_n, mirror_point

logger = logging.getLogger('femagtools.symmetry')


#############################
#         symmetry          #
#############################


class Symmetry(object):
    def __init__(self,
                 geom=None,
                 startangle=None,
                 endangle=None,
                 rtol=1e-04,
                 atol=1e-03):
        assert(geom is not None)
        assert(startangle is not None)
        assert(endangle is not None)
        
        self.geom = geom
        self.startangle = startangle
        self.endangle = endangle
        self.delta_check_count = 0
        self.delta_angle_korr = 0.0
        self.rtol = rtol
        self.atol = atol
        self.full = False
        if np.isclose(self.startangle, self.endangle):
            self.alpha = 2.0*np.pi
            self.full = True
        else:
            self.alpha = alpha_angle(self.startangle, self.endangle)
            self.full = False
        logger.debug("Symmetry(alpha=%s, rtol=%s, atol=%s)", self.alpha, rtol, atol)

    def __str__(self):
        return "rtol: {}\n".format(self.rtol) + \
               "atol: {}\n".format(self.atol)

    def equal_area(self,
                   d1, h1, a1,
                   d2, h2, a2,
                   rtol=1e-03, atol=1e-03):
        if not np.isclose(d1, d2, rtol=rtol, atol=atol):
            logger.debug("dist NOT close (%s/%s)", d1, d2)
            return False
        if not np.isclose(h1, h2, rtol=rtol, atol=atol):
            logger.debug("height NOT close (%s/%s)", h1, h2)
            return False
        if not np.isclose(a1, a2, rtol=rtol, atol=atol):
            logger.debug("alpha NOT close (%s/%s)", a1, a2)
            return False
        return True

    def calc_mid_angle(self, a):
        return positive_angle(alpha_angle(self.startangle,
                                          a.get_mid_angle(self.geom.center)))

    def find_symmetry(self):
        arealist = self.geom.list_of_areas()
        logger.debug("begin of Symmetry::find_symmetry: %s areas available", len(arealist))
        if len(arealist) == 0:
            logger.debug("end of find_symmetry: no areas")
            return 0

        areas = []
        for a in arealist:
            areas.append((round(a.get_alpha(self.geom.center), 3),
                          round(a.min_dist, 1),
                          round(a.height, 1),
                          self.calc_mid_angle(a),
                          a))
        areas.sort(reverse=True)

        a0_alpha, a0_min_dist, a0_height, a0_mid_angle, a0 = areas[0]
        equal_areas = [(a0_mid_angle, a0)]
        check_rslt = []
        for a1_alpha, a1_min_dist, a1_height, a1_mid_angle, a1 in areas[1:]:
            if self.equal_area(a0_min_dist, a0_height, a0_alpha,
                               a1_min_dist, a1_height, a1_alpha,
                               rtol=0.01, atol=0.05):
                a0_min_dist = (a0_min_dist + a1_min_dist) / 2
                a0_height = (a0_height + a1_height) / 2
                a0_alpha = (a0_alpha + a1_alpha) / 2
                equal_areas.append((a1_mid_angle, a1))
            else:
                parts = self.check_delta(equal_areas)
                check_rslt.append((a0.area_size(), parts, len(equal_areas), a0))

                equal_areas = [(a1_mid_angle, a1)]
                a0_min_dist = a1_min_dist
                a0_height = a1_height
                a0_alpha = a1_alpha
                a0 = a1

        parts = self.check_delta(equal_areas)
        check_rslt.append((a0.area_size(), parts, len(equal_areas), a0))

        parts = self.get_symmetry_parts(check_rslt)
        if parts < 2:
            logger.debug("end of Symmetry::find_symmetry: no symmetry")
            return parts

        self.geom.clear_cut_lines()
        for alpha in self.symmetry_lines(parts,
                                         self.startangle,
                                         self.endangle):
            plus = self.geom.max_radius / 10
            min_radius = max(10, self.geom.min_radius - plus)
            p1 = point(self.geom.center, min_radius, alpha)
            p2 = point(self.geom.center, self.geom.max_radius + plus, alpha)
            line = Line(Element(start=p1, end=p2))
            line.init_attributes(color='green')
            self.geom.add_cut_line(line)

        logger.debug("end of Symmetry::find_symmetry: -> %s", parts)
        return parts
    
    def check_delta(self, area_list):
        logger.debug("begin of check_delta: %s equal areas", len(area_list))
        if not area_list:
            logger.debug("end of check_delta: no areas")
            return None

        if len(area_list) == 1:
            mid_angle, a = area_list[0]
            alpha = a.get_alpha(self.geom.center)
            if np.isclose(alpha, self.alpha):
                logger.debug("end of check_delta: one area from start to end")
                return None  # ok

        self.delta_check_count += 1
        area_list.sort()

        mid_angle, a = area_list[0]
        delta = positive_angle(mid_angle * 2)
        delta_total = mid_angle
        delta_list = [delta]

        logger.debug("First mid = %s, delta = %s", mid_angle, delta)
        logger.debug("%s:  d=%s,  h=%s,  a=%s, mid=%s, delta=%s",
                     a.identifier(),
                     a.min_dist,
                     a.height,
                     a.get_alpha(self.geom.center),
                     mid_angle,
                     delta)

        geom_alpha = alpha_angle(self.startangle, self.endangle)
        geom_alpha = positive_angle(geom_alpha)

        start_angle = mid_angle
        for mid_angle, a in area_list[1:]:
            delta_angle = alpha_angle(start_angle, mid_angle)
            delta = positive_angle(delta_angle)
            delta_total += delta_angle
            delta_list.append(delta)

            logger.debug("%s:  d=%s,  h=%s,  a=%s, mid=%s, delta=%s",
                         a.identifier(),
                         a.min_dist,
                         a.height,
                         a.get_alpha(self.geom.center),
                         mid_angle,
                         delta)
            start_angle = mid_angle
            
        delta_angle = alpha_angle(start_angle, geom_alpha)
        delta = positive_angle(delta_angle * 2)
        delta_total += delta_angle
        delta_list.append(delta)
        logger.debug("final delta=%s", delta)

        if not np.isclose(geom_alpha, delta_total):
            logger.debug("-- deltas: %s", delta_list)
            logger.debug("end of check_delta: BAD DELTA %s, (expected %s)",
                         delta_angle, geom_alpha)
            return 0  # very bad

        sz = len(delta_list)
        mid = int(sz / 2)
        ix1 = 0
        ix2 = sz - 1
        first_last_bad = False
        for ix1 in range(0, mid):
            if not np.isclose(delta_list[ix1], delta_list[ix2], rtol=1e-3, atol=1e-3):
                if self.full and \
                   self.delta_check_count == 1 and \
                   ix1 == 0 and \
                   self.delta_angle_korr == 0.0:
                    first_last_bad = True
                else:
                    logger.debug("end of check_delta: NO SYM")
                    return 0
            ix2 -= 1

        if first_last_bad:
            delta_korr = (delta_list[0] + delta_list[-1]) / 2.0
            logger.debug("STARTANGLE CORRECTION")
            self.delta_angle_korr = (delta_korr - delta_list[0]) / 2
            logger.debug("-- delta[0] from %s to %s", delta_list[0], delta_korr)
            logger.debug("Delta Angle Korr = %s", self.delta_angle_korr)
            delta_list[0] = delta_korr
            delta_list[-1] = delta_korr
            assert(self.full)
            self.startangle = self.startangle - self.delta_angle_korr
            self.endangle = self.endangle - self.delta_angle_korr
            logger.debug("New startangle = %s", self.startangle)
            logger.debug("Delta List: %s", delta_list)

        d1 = delta_list[0]
        d1_count = 1
        inx_list = [0]
        for x in range(1, len(delta_list)):
            if np.isclose(d1, delta_list[x], rtol=1e-3, atol=1e-3):
                inx_list.append(x)
                d1_count += 1

        if d1_count == len(delta_list):
            logger.debug("end of check_delta: SYMMETRY FOUND")
            return d1_count -1  # very simple
        if len(delta_list) < 2:
            logger.debug("end of check_delta: One delta only ?!")
            return 0

        logger.debug("index of delta %s = %s", d1, inx_list)
        x1 = inx_list[0]
        x2 = inx_list[1]
        step = x2 - x1
        x1 = x2
        for x2 in inx_list[2:]:
            if not (x2 - x1 == step):
                return 0
            x1 = x2
        
        logger.debug("end of check_delta: SYMMETRY FOUND")
        return len(inx_list) -1

    def get_symmetry_parts(self, check_rslt):
        max_size = 0
        max_areas = 0
        parts_possible = None

        check_rslt.sort(reverse=True)
        for size, parts, count, area in check_rslt:
            logger.debug("Result: %s, %s, %s", size, parts, count)

        for size, parts, count, area in check_rslt:
            if parts is not None and parts > 0:
                max_size = max(max_size, size)
                max_areas = max(max_areas, count)

        logger.debug("max size: %s,  max areas: %s", max_size, max_areas)

        for size, parts, count, area in check_rslt:
            if parts is not None and parts <= 1:  # critical
                if count <= max(1, max_areas / 5):
                    if size < max_size / 25:
                        parts = None

            parts_possible = self.calc_parts(parts_possible, parts)

        if parts_possible is None:
            parts_possible = 0
        return parts_possible

    def calc_parts(self, parts1, parts2):
        logger.debug("Calc symmetry Parts (%s, %s)", parts1, parts2)
        if parts2 is None:
            logger.debug("return %s parts", parts1)
            return parts1
        if parts1 is None:
            logger.debug("return %s parts", parts2)
            return parts2
        if parts1 == 0 or parts2 == 0:
            logger.debug("return %s parts", 0)
            return 0
        parts = gcd(parts1, parts2)
        logger.debug("return %s parts", parts)
        return parts

    def symmetry_lines(self, parts, startangle, endangle):
        logger.debug("begin symmetry_lines from %s to %s",
                     startangle,
                     endangle)
        if less_equal(endangle, startangle):
            endangle += 2*np.pi

        delta = alpha_angle(startangle, endangle) / parts
        start = startangle + delta
        while less(start, endangle):
            yield start
            start += delta

        if is_same_angle(startangle, endangle):
            yield start

        # Damit man anschliessend ohne Umstände schneiden kann.
        self.geom.sym_startangle = startangle
        self.geom.sym_endangle = startangle + delta
        self.geom.sym_slices = parts
        self.geom.sym_slice_angle = delta
        self.geom.sym_area = Area([], (0,0), 0.0)
        self.geom.sym_area.sym_startangle = self.geom.sym_startangle
        self.geom.sym_area.sym_endangle = self.geom.sym_endangle
        logger.debug("end symmetry_lines")

    def check_symmetry_of_mirror(self, mirror_geom, mirrorangle):
        logger.debug("begin of Symmetry::check_symmetry_of_mirror")
        assert(mirror_geom is not None)

        axis_p = point(self.geom.center, self.geom.max_radius, mirrorangle)
        axis_m = line_m(self.geom.center, axis_p)
        axis_n = line_n(self.geom.center, axis_m)

        def counterpart_found(node, nodes, rtol, atol):
            hits = 0
            for n in nodes:
                if points_are_close(node, n, rtol, atol):
                    logger.debug(" ---- %s is %s", node, n)
                    return True
            return False

        def check_differences(geom, mirror_geom):
            geom_ag_nodes = []
            geom_nodes = [n for n in geom.g.nodes() if not (n in geom_ag_nodes)]

            hit = 0
            for n in geom_nodes:
                mirror_n = mirror_point(n, geom.center, axis_m, axis_n)
                if counterpart_found(mirror_n,
                                     mirror_geom.g.nodes(),
                                     self.rtol,
                                     self.atol):
                    hit += 1
            min_nodes = min(len(geom_nodes), int(len(geom_nodes) * 0.95) + 1)
            logger.debug("Nodes=%s,  Counterparts=%s", len(geom_nodes), hit)
            if hit < min_nodes:
                return hit / len(geom_nodes)
            else:
                return 1.0

        # ----------------
        logger.debug("check geom - mirror")
        f1 = check_differences(self.geom, mirror_geom)
        logger.debug("check mirror - geom")
        f2 = check_differences(mirror_geom, self.geom)
        logger.debug("Factor 1: %s,  2: %s", f1, f2)
        ok = f1 > 0.97 and f2 > 0.97
        logger.debug("end of Symmetry::check_symmetry_of_mirror => %s", ok)
        return ok
