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
from femagtools.dxfsl.functions import alpha_angle, positive_angle
from femagtools.dxfsl.functions import min_angle, max_angle, gcd, point
from femagtools.dxfsl.functions import less_equal, less

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
        self.rtol = rtol
        self.atol = atol

        logger.debug("Symmetry(rtol=%s, atol=%s)", rtol, atol)

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

    def find_symmetry(self):
        arealist = self.geom.list_of_areas()
        logger.debug("begin of Symmetry::find_symmetry: %s areas available", len(arealist))
        if len(arealist) == 0:
            logger.debug("end of find_symmetry: no areas")
            return 0

        areas = []
        for a in arealist:
            areas.append((a.get_alpha(),
                          a.get_mid_angle(),
                          a.min_dist,
                          a.height,
                          a))
        areas.sort()
        parts_possible = None

        a0_alpha, a0_mid_angle, a0_min_dist, a0_height, a0 = areas[0]
        equal_areas = [(a0_mid_angle, a0)]

        for a1_alpha, a1_mid_angle, a1_min_dist, a1_height, a1 in areas[1:]:
            if self.equal_area(a0_min_dist, a0_height, a0_alpha,
                               a1_min_dist, a1_height, a1_alpha,
                               rtol=0.01, atol=0.05):
                a0_min_dist = (a0_min_dist + a1_min_dist) / 2
                a0_height = (a0_height + a1_height) / 2
                a0_alpha = (a0_alpha + a1_alpha) / 2
                equal_areas.append((a1_mid_angle, a1))
            else:
                parts = self.check_delta(equal_areas)
                if parts_possible is None:
                    parts_possible = parts
                else:
                    parts_possible = self.calc_parts(parts_possible, parts)
                equal_areas = [(a1_mid_angle, a1)]
                a0_min_dist = a1_min_dist
                a0_height = a1_height
                a0_alpha = a1_alpha
                a0 = a1

        parts = self.check_delta(equal_areas)
        if parts_possible is None:
            parts_possible = parts
        else:
            parts_possible = self.calc_parts(parts_possible, parts)

        if parts_possible is None:
            parts_possible = 0
        if parts_possible < 2:
            logger.debug("end of Symmetry::find_symmetry: no symmetry")
            return parts_possible
    
        self.geom.clear_cut_lines()
        for alpha in self.symmetry_lines(parts_possible,
                                         self.startangle,
                                         self.endangle):
            min_radius = max(10, self.geom.min_radius - 10)
            p1 = point(self.geom.center, min_radius, alpha)
            p2 = point(self.geom.center, self.geom.max_radius+10, alpha)
            line = Line(Element(start=p1, end=p2))
            line.init_attributes(color='green')
            self.geom.add_cut_line(line)
                    
        logger.debug("end of Symmetry::find_symmetry")
        return parts_possible
    
    def check_delta(self, area_list):
        logger.debug("==> %s equal areas", len(area_list))
        if not area_list:
            return None

        if len(area_list) == 1:
            mid_angle, a = area_list[0]
            if np.isclose(a.min_angle, self.startangle) and \
               np.isclose(a.max_angle, self.endangle):
                return None  # ok

        area_list.sort()
        start = self.startangle
        mid_angle, a = area_list[0]

        delta = alpha_angle(start, mid_angle) * 2
        delta_list = [delta]
        logger.debug("geom: start=%s,  end=%s", self.startangle, self.endangle)
        logger.debug("%s:  d=%s,  h=%s,  a=%s, mid=%s, delta=%s",
                     a.identifier(),
                     a.min_dist,
                     a.height,
                     a.get_alpha(),
                     a.get_mid_angle(),
                     delta)
        logger.debug("  min=%s,  max%s",
                     a.min_angle,
                     a.max_angle)

        start = mid_angle
        for mid_angle, a in area_list[1:]:
            delta = alpha_angle(start, mid_angle)
            delta_list.append(delta)

            logger.debug("%s:  d=%s,  h=%s,  a=%s, mid=%s, delta=%s",
                         a.identifier(),
                         a.min_dist,
                         a.height,
                         a.get_alpha(),
                         a.get_mid_angle(),
                         delta)
            logger.debug("  min=%s,  max%s",
                         a.min_angle,
                         a.max_angle)

            start = mid_angle
            
        delta = alpha_angle(start, self.endangle) * 2
        delta_list.append(delta)
        logger.debug("final delta=%s", delta)

        sz = len(delta_list)
        mid = int(sz / 2)
        ix1 = 0
        ix2 = sz - 1
        for ix1 in range(0, mid):
            if not np.isclose(delta_list[ix1], delta_list[ix2], rtol=1e-3, atol=1e-3):
                logger.debug("NO SYM")
                return 0
            ix2 -= 1

        geom_alpha = alpha_angle(self.startangle, self.endangle)
        delta_alpha = 0.0
        for d in delta_list:
            delta_alpha += d

        geom_alpha = positive_angle(geom_alpha)

        delta_alpha -= delta_list[0] / 2
        delta_alpha -= delta_list[-1] / 2
        delta_alpha = positive_angle(delta_alpha)
        if not np.isclose(geom_alpha, delta_alpha):
            logger.debug("BAD DELTA ?!?")
            return 0  # very bad

        d1 = delta_list[0]
        d1_count = 1
        inx_list = [0]
        for x in range(1, len(delta_list)):
            if np.isclose(d1, delta_list[x], rtol=1e-3, atol=1e-3):
                inx_list.append(x)
                d1_count += 1

        if d1_count == len(delta_list):
            logger.debug("SYMMETRY FOUND")
            return d1_count -1  # very simple
        if len(delta_list) < 2:
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
        
        logger.debug("SYMMETRY FOUND")
        return len(inx_list) -1

    def calc_parts(self, parts1, parts2):
        logger.debug("Calc symmetry Parts (%s, %s)", parts1, parts2)
        if parts2 is None:
            return parts1

        return gcd(parts1, parts2)

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
            
        # Damit man anschliessend ohne Umstände schneiden kann.
        self.geom.sym_startangle = startangle
        self.geom.sym_endangle = startangle + delta
        self.geom.sym_slices = parts
        self.geom.sym_slice_angle = delta
        self.geom.sym_area = Area([], (0,0), 0.0)
        self.geom.sym_area.sym_startangle = self.geom.sym_startangle
        self.geom.sym_area.sym_endangle = self.geom.sym_endangle
        logger.debug("end symmetry_lines")
