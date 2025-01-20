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
import femagtools.dxfsl.area as AREA
from femagtools.dxfsl.functions import alpha_angle, positive_angle, is_same_angle
from femagtools.dxfsl.functions import min_angle, max_angle, gcd, point
from femagtools.dxfsl.functions import less_equal, less, points_are_close
from femagtools.dxfsl.functions import line_m, line_n, mirror_point
from femagtools.dxfsl.functions import part_of_circle, round_point

logger = logging.getLogger('femagtools.symmetry')


#############################
#         symmetry          #
#############################


class Symmetry(object):
    def __init__(self,
                 geom=None,
                 startangle=None,
                 endangle=None,
                 rtol=1e-03,
                 atol=1e-02):
        assert(geom is not None)
        assert(startangle is not None)
        assert(endangle is not None)
        
        self.geom = geom
        self.startangle = startangle
        self.endangle = endangle
        self.geom_part = part_of_circle(self.startangle, self.endangle, 1)
        self.delta_check_count = 0
        self.delta_angle_corr = None
        self.rtol = rtol
        self.atol = atol
        self.full = False
        self.ag_radius = 0.0
        if np.isclose(self.startangle, self.endangle):
            self.alpha = 2.0*np.pi
            self.full = True
        else:
            self.alpha = alpha_angle(self.startangle, self.endangle)
            self.full = False
        if geom.is_inner:
            self.ag_radius = geom.max_radius
        else:
            self.ag_radius = geom.min_radius
        logger.debug("Symmetry(alpha=%s, rtol=%s, atol=%s)", self.alpha, rtol, atol)

    def __str__(self):
        return "rtol: {}\n".format(self.rtol) + \
               "atol: {}\n".format(self.atol)

    def equal_area(self,
                   d1, h1, a1,
                   d2, h2, a2,
                   rtol=1e-03, atol=1e-03):
        equal_d = np.isclose(d1, d2, rtol=rtol, atol=atol)  # distance form c
        equal_h = np.isclose(h1, h2, rtol=rtol, atol=atol)  # height
        equal_a = np.isclose(a1, a2, rtol=1e-4, atol=1e-2)  # angle from c
        if not equal_d:
            logger.debug("equal_area: dist NOT close (%s/%s)", d1, d2)
            if equal_h and equal_a:  # but height and angle
                if np.isclose(d1, d2, rtol=1e-2, atol=1.e-1):
                    logger.debug(" -- but with more tolerance")
                    return True
            return False
        if not equal_h:
            logger.debug("equal_area: height NOT close (%s/%s)", h1, h2)
            if equal_d and equal_a:  # but distance and angle
                if np.isclose(h1, h2, rtol=1e-2, atol=1e-1) :
                    logger.debug(" -- but with more tolerance")
                    return True
            return False
        if not np.isclose(a1, a2, rtol=1e-4, atol=1e-2):
            logger.debug("equal_area: alpha NOT close (%s/%s)", a1, a2)
            return False
        else:
            if a1 > a2:
                f = a2 / a1
            else:
                f = a1 / a2
            if f < 0.9:
                return False
        return True

    def calc_mid_angle(self, a):
        return positive_angle(alpha_angle(self.startangle,
                                          a.get_mid_angle(self.geom.center)))

    def area_list_entry(self, a):
        a.set_symmetry_parameter(self.geom.center)
        return (round(a.get_alpha(self.geom.center), 3),
                round(a.min_dist, 1),
                round(a.height, 1),
                self.calc_mid_angle(a),
                a)

    def build_area_list(self, types=()):
        arealist = self.geom.list_of_areas()
        if types:
            arealist = [a for a in arealist if a.type in types]

        areas = []
        for a in arealist:
            areas.append(self.area_list_entry(a))
        areas.sort(reverse=True)
        return areas

    def get_equal_areas(self, areas):
        return

    def build_results(self, areas):
        logger.debug("begin of build_results with %s areas", len(areas))
        [logger.debug("#%s: alpha=%s, min=%s, max=%s",
                      a.get_id(), alpha, a.min_angle, a.max_angle) for
         alpha, dist, h, mid, a in areas]

        a0_alpha, a0_min_dist, a0_height, a0_mid_angle, a0 = areas[0]
        equal_areas = [(a0_mid_angle, a0)]
        check_rslt = []
        for a1_alpha, a1_min_dist, a1_height, a1_mid_angle, a1 in areas[1:]:
            if (self.equal_area(a0_min_dist, a0_height, a0_alpha,
                                a1_min_dist, a1_height, a1_alpha,
                                rtol=0.001, atol=0.05)):
                a0_min_dist = (a0_min_dist + a1_min_dist) / 2
                a0_height = (a0_height + a1_height) / 2
                a0_alpha = (a0_alpha + a1_alpha) / 2
                equal_areas.append((a1_mid_angle, a1))
            else:
                # alpha Wechsel
                id_list = []
                for i in range(len(equal_areas)):
                    mid_angle0, area0 = equal_areas[i]
                    if area0.get_id() in id_list:
                        continue
                    equal_areas_check = [(mid_angle0, area0)]
                    for mid_angle1, area1 in equal_areas[i+1:]:
                        if area1.get_id() in id_list:
                            continue
                        if area0.is_symmetry_equal(area1):
                            equal_areas_check.append((mid_angle1, area1))
                            id_list.append(area1.get_id())

                    rslt = self.check_delta(equal_areas_check)
                    areasize = a0.area_size()
                    rslt['area'] = a0
                    rslt['areasize'] = areasize
                    check_rslt.append((areasize, rslt))

                equal_areas = [(a1_mid_angle, a1)]
                a0_min_dist = a1_min_dist
                a0_height = a1_height
                a0_alpha = a1_alpha
                a0 = a1

        rslt = self.check_delta(equal_areas)
        areasize = a0.area_size()
        rslt['area'] = a0
        rslt['areasize'] = areasize
        check_rslt.append((areasize, rslt))
        logger.debug("end of build_results")
        return check_rslt

    def get_winding_symmetry(self, inside=False):
        if inside:
            areas = [self.area_list_entry(a) for a in self.geom.list_of_areas()
                     if not a.close_to_ag]
            areas.sort(reverse=True)
        else:
            areas = self.build_area_list((AREA.TYPE_WINDINGS,))

        logger.debug("begin of Symmetry::get_winding_symmetry: %s areas available", len(areas))
        if not areas:
            logger.debug("end of Symmetry::get_winding_symmetry: no areas")
            return 0

        check_rslt = self.build_results(areas)
        logger.debug("%s results available", len(check_rslt))
        [logger.debug("Result: %s", rslt) for rslt in check_rslt]

        parts, start_delta = self.get_symmetry_parts(check_rslt)
        if parts <= 1:
            return 0
        self.create_cut_lines(parts, start_delta)

        sym = self.geom_part * parts
        delta = 2*np.pi/sym
        self.set_symmetry_parameters(self.startangle, parts, delta)

        logger.debug("end of Symmetry::get_winding_symmetry: parts=%s", parts)
        return parts

    def get_magnet_symmetry(self):
        areas = self.build_area_list((AREA.TYPE_MAGNET_AIRGAP, AREA.TYPE_MAGNET_RECT,))
        air = self.build_area_list((AREA.TYPE_AIR,))
        mag_list = [a for a in self.geom.list_of_areas() if a.is_magnet()]
        air_list = [a for a in self.geom.list_of_areas() if a.is_air()]
        sz_list = [a.area_size() for a in self.geom.list_of_areas()]
        max_sz = max(sz_list)
        for a in air_list:
            if a.area_size() < max_sz * 0.005:
                continue
            for m in mag_list:
                if a.is_touching(m):
                    areas.append(self.area_list_entry(a))
                    break

        logger.debug("begin of Symmetry::get_magnet_symmetry: %s areas available", len(areas))
        if not areas:
            logger.debug("end of Symmetry::get_magnet_symmetry: no areas")
            return 0

        check_rslt = self.build_results(areas)
        logger.debug("%s results available", len(check_rslt))
        [logger.debug("Result: %s", rslt) for rslt in check_rslt]
        for sz, rslt in check_rslt:
            if not rslt.get('startdelta', 0.0) == 0.0:
                return 0  # not proper
            if rslt.get('halfslice', None):
                return 0  # not proper

        parts, start_delta = self.get_symmetry_parts(check_rslt)
        if parts <= 1:
            return 0
        self.create_cut_lines(parts, start_delta)

        sym = self.geom_part * parts
        delta = 2*np.pi/sym
        self.set_symmetry_parameters(self.startangle, parts, delta)

        logger.debug("end of Symmetry::get_magnet_symmetry: parts=%s", parts)
        return parts

    def find_symmetry(self):
        areas = self.build_area_list()

        logger.debug("begin of Symmetry::find_symmetry: %s areas available", len(areas))
        if not areas:
            logger.debug("end of Symmetry::find_symmetry: no areas")
            return 0

        check_rslt = self.build_results(areas)
        logger.debug("%s results available", len(check_rslt))

        parts, start_delta = self.get_symmetry_parts(check_rslt)
        if parts < 2:
            logger.debug("end of Symmetry::find_symmetry: no symmetry")
            return parts

        if self.delta_angle_corr is not None and self.delta_angle_corr != 0.0:
            self.startangle = self.startangle - self.delta_angle_corr
            self.endangle = self.endangle - self.delta_angle_corr

        self.create_cut_lines(parts, start_delta)

        logger.debug("end of Symmetry::find_symmetry: -> %s", parts)
        return parts

    def check_delta(self, area_list):
        logger.debug("begin of check_delta: %s equal areas", len(area_list))
        result = {'areas': len(area_list),
                  'startdelta': 0.0,
                  'slices': None}
        result['area_id_list'] = [a.get_id() for m, a in area_list]
        if not area_list:
            logger.debug("end of check_delta: no areas")
            return result

        rtol = 1e-3
        atol = 1e-2

        logger.debug("Geometry: Alpha=%s,  Center=%s", self.alpha, self.geom.center)
        mid_angle, a = area_list[0]
        result['height'] = a.height
        result['alpha'] = a.get_alpha(self.geom.center)
        if self.geom.is_inner:
            result['airgap'] = np.isclose(a.max_dist, self.ag_radius)
        else:
            result['airgap'] = np.isclose(a.min_dist, self.ag_radius)
        if len(area_list) == 1:  # one area
            return self.check_one_area(mid_angle, a, result, rtol=rtol, atol=atol)

        self.delta_check_count += 1
        area_list.sort()

        mid_angle, a = area_list[0]
        delta = positive_angle(mid_angle * 2)
        delta_total = mid_angle
        logger.debug("first delta=%s,  total=%s", delta, delta_total)
        delta_list = [delta]
        mid_delta_list = [(mid_angle, delta)]
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
            logger.debug("next delta=%s,  total=%s", delta, delta_total)

            delta_list.append(delta)
            mid_delta_list.append((mid_angle, delta))

            logger.debug("%s:  d=%s,  h=%s,  a=%s, mid=%s, delta=%s",
                         a.identifier(),
                         a.min_dist,
                         a.height,
                         a.get_alpha(self.geom.center),
                         mid_angle,
                         delta)
            start_angle = mid_angle
            
        delta_angle = alpha_angle(start_angle, geom_alpha)
        if np.isclose(delta_angle, np.pi*2, rtol=1e-4, atol=1e-4):
            logger.debug("Last Area is in the middle of endangle")
            delta_angle = 0.0
        delta = positive_angle(delta_angle * 2)
        delta_total += delta_angle
        delta_list.append(delta)
        mid_delta_list.append((0.0, delta))

        logger.debug("final delta=%s,  total=%s", delta, delta_total)
        logger.debug("Delta List: %s", delta_list)
        logger.debug("Mid Delta List")
        [logger.debug("-- Mid angle: %s,   Delta: %s", a, d) for a, d in mid_delta_list]

        if not np.isclose(geom_alpha, delta_total, rtol=rtol, atol=atol):
            logger.debug("-- deltas: %s", delta_list)
            logger.debug("end of check_delta: BAD DELTA %s, (expected %s)",
                         delta_angle, geom_alpha)
            result['slices'] = 0
            return result  # very bad

        deltas = self.create_deltas(delta_list, rtol=rtol, atol=atol)

        logger.debug("Start with looking for symmetry")

        logger.debug(">> %s Deltas <<", len(deltas))
        [logger.debug(" -- n=%s,  delta=%s", n, d) for n, d in deltas]

        if len(deltas) == 2:
            n1, d1 = deltas[0]
            n2, d2 = deltas[1]
            logger.debug("delta 1: n=%s,  delta=%s", n1, d1)
            logger.debug("delta 2: n=%s,  delta=%s", n2, d2)

            if n2 == 2 and n1 > 2 and \
               np.isclose(d1, d2 / 2.0, rtol=rtol, atol=atol):
                if np.isclose(d2, delta_list[0], rtol=rtol, atol=atol) and \
                   np.isclose(d2, delta_list[-1], rtol=rtol, atol=atol):
                    slices = n1 + n2
                    if slices > 4:
                        result['slices'] = slices
                    else:
                        result['slices'] = 1
                    result['slices_half'] = slices
                    result['halfslice'] = 1
                    result['startdelta'] = d1 / 2
                    logger.debug("#3: end of check_delta: SYMMETRY FOUND [%s] halfslice",
                                 result['slices'])
                    return result

            elif n1 == 2 and n2 == 1:
                assert(len(area_list) == 2)
                semi_alpha = result['alpha'] / 2
                a0_mid_angle, a0 = area_list[0]
                a1_mid_angle, a1 = area_list[1]
                a0_start = np.isclose(a0_mid_angle - semi_alpha,
                                      0.0,
                                      rtol=rtol, atol=atol)
                a1_end = np.isclose(a1_mid_angle + semi_alpha,
                                    self.alpha,
                                    rtol=rtol, atol=atol)
                if a0_start and a1_end and d1 < d2:
                    result['slices'] = 0
                    result['halfslice'] = 2
                    logger.debug("#4: end of check_delta: half slices")
                    return result

                if not a0_start and not a1_end and d1 > d2:
                    parts_in_circ = part_of_circle(0.0, d2, 1)
                    parts_in_geom = float(round(parts_in_circ / self.geom_part, 2))
                    if parts_in_geom.is_integer():
                        parts_in_geom = int(parts_in_geom)
                    else:
                        parts_in_geom = 0
                    result['slices'] = parts_in_geom
                    result['slices_half'] = parts_in_geom
                    result['halfslice'] = 1
                    logger.debug("#5: end of check_delta: SYMMETRY FOUND [%s] halfslice",
                                 result['slices'])
                    return result

            elif abs(n1 - n2) == 1:
                if d1 < d2:
                    if self.check_pairs_in_delta_list(delta_list,
                                                      deltas,
                                                      rtol=rtol, atol=atol):
                        result['slices'] = int(len(area_list) / 2)
                        delta_angle_corr = (d1 + d2) / 2
                        result['delta_corr'] = delta_angle_corr
                        logger.debug("Startangle correction by %s", delta_angle_corr)
                        logger.debug("#6: end of check_delta: SYMMETRY FOUND [%s]",
                                     result['slices'])
                        return result

        if len(deltas) > 1:
            n1, d1 = deltas[0]
            n2, d2 = deltas[1]
            if n1 + n2 + 1 == len(area_list) and abs(n1 - n2) == 1:
                dlist = self.check_pairs_of_areas(mid_delta_list,
                                                  min(d1, d2),
                                                  geom_alpha,
                                                  rtol=rtol, atol=atol)
                if dlist:
                    delta_list = dlist
                    deltas = self.create_deltas(delta_list, rtol=rtol, atol=atol)

        if False:
            parts_in_circ = part_of_circle(0.0, d1, 1)
            parts_in_geom = float(round(parts_in_circ / self.geom_part, 2))
            if parts_in_geom.is_integer():
                parts_in_geom = int(parts_in_geom)
            else:
                parts_in_geom = 0

            if parts_in_geom / n1 > 0.75 and parts_in_circ > 15:
                result['slices'] = parts_in_geom
                logger.debug("#7: end of check_delta: SYMMETRY FOUND [%s]",
                             result['slices'])
                return result

        missing_middle = []
        if len(deltas) == 2:
            logger.debug("looking for holes in delta list")

            delta_n, delta_value = deltas[0]
            logger.debug("First n,v == (%s, %s)", delta_n, delta_value)
            for n, v in deltas[1:]:
                logger.debug("Next n,v == (%s, %s)", n, v)
                if n < delta_n / 4 and \
                   np.isclose(delta_value, v / 2.0, rtol=1e-04, atol=1e-03):
                    logger.debug("Hole found")
                    inx = [i for i, x in enumerate(delta_list)
                           if np.isclose(x, v, rtol=rtol, atol=atol)]
                    if len(inx) != n:
                        logger.debug("Hole missmatch: %s <> %s", len(inx), n)
                        result['slices'] = 0
                        return result

                    dlist = []
                    x = 0

                    for i in inx:
                        for n in range(x, i):
                            logger.debug("set value of index %s", n)
                            dlist.append(delta_list[n])
                        m1, a = area_list[i-1]
                        if i < len(area_list):
                            m2, a = area_list[i]
                        else:
                            m2, a = area_list[0]
                        mid = (m1 + m2) / 2
                        logger.debug("Missing mid is %s", mid)
                        missing_middle.append((mid, i))
                        x = i+1
                        logger.debug("set value in hole")
                        dlist.append(delta_value)
                        dlist.append(delta_value)
                    for n in range(x, len(delta_list)):
                        logger.debug("set value of index %s", n)
                        dlist.append(delta_list[n])
                    logger.debug("New List: %s", dlist)
                    delta_list = dlist

        result['middlelist'] = [m for m, a in area_list]
        result['missing_middles'] = missing_middle

        if not np.isclose(delta_list[0], delta_list[-1], rtol=rtol, atol=atol):
            logger.debug("First and Last delta not equal")
            if self.full:
                d0 = (delta_list[0] + delta_list[-1]) / 2
                n1, d1 = deltas[0]
                if np.isclose(d0, d1, rtol=rtol, atol=atol):
                    delta_angle_corr = (d0 - delta_list[0]) / 2
                    result['delta_corr'] = delta_angle_corr
                    logger.debug("Startangle correction by %s", delta_angle_corr)
                    delta_list[0] = d0
                    delta_list[-1] = d0
            else:
                parts = self.check_first_last_difference(delta_list, deltas)
                if parts > 0:
                    result['slices'] = parts
                    logger.debug("#8: end of check_delta: SYMMETRY FOUND [%s]",
                                 result['slices'])
                    return result

        logger.debug("Final Delta List: %s", delta_list)
        d1 = delta_list[0]
        d1_count = 1
        inx_list = [0]
        for x in range(1, len(delta_list)):
            if np.isclose(d1, delta_list[x], rtol=rtol, atol=atol):
                inx_list.append(x)
                d1_count += 1

        if d1_count == len(delta_list):
            result['slices'] = d1_count -1
            logger.debug("#9: end of check_delta: SYMMETRY FOUND [%s]",
                         result['slices'])
            return result  # very simple
        if len(delta_list) < 3:
            logger.debug("end of check_delta: One delta only ?!")
            result['slices'] = 0
            return result

        logger.debug("index of delta %s: %s", d1, inx_list)
        if len(inx_list) < 2:
            logger.debug("end of check_delta: NO SYMMETRY")
            result['slices'] = 0
            return result

        x1 = inx_list[0]
        x2 = inx_list[1]
        step = x2 - x1
        x1 = x2
        for x2 in inx_list[2:]:
            if not (x2 - x1 == step):
                logger.debug("end of check_delta: NO SYMMETRY")
                result['slices'] = 0
                return result
            x1 = x2

        logger.debug("length of delta %s: %s",
                     len(inx_list),
                     inx_list)
        result['slices'] = len(inx_list) -1
        logger.debug("#10: end of check_delta: SYMMETRY FOUND [%s]", result['slices'])
        return result

    def create_deltas(self, delta_list, rtol=1e-3, atol=1e-2):
        delta_list_sorted = [d for d in delta_list]
        delta_list_sorted.sort()
        delta = delta_list_sorted[0]
        delta_n = 1
        deltas = []
        for i in range(1,len(delta_list_sorted)):
            if np.isclose(delta_list_sorted[i], delta, rtol=rtol, atol=atol):
                delta = (delta + delta_list_sorted[i]) / 2
                delta_n += 1
            else:
                deltas.append((delta_n, delta))
                delta = delta_list_sorted[i]
                delta_n = 1
        deltas.append((delta_n, delta))
        deltas.sort(reverse=True)
        return deltas

    def check_one_area(self, mid_angle, a, result, rtol=1e-3, atol=1e-2):
        logger.debug("begin of check_one_area")

        alpha = a.get_alpha(self.geom.center)
        logger.debug("Single %s:  d=%s,  h=%s,  a=%s, mid=%s",
                     a.identifier(),
                     a.min_dist,
                     a.height,
                     a.get_alpha(self.geom.center),
                     mid_angle)
        if np.isclose(alpha, self.alpha, rtol=rtol, atol=atol):
            logger.debug("end of check_one_area: area %s from start to end",
                         a.identifier())
            result['slices'] = None
            return result  # ok

        if self.full:
            result['slices'] = 1
            logger.debug("end of check_one_area: full with 1 slice")
            return result

        delta_angle = alpha_angle(self.startangle, mid_angle)
        delta1 = positive_angle(delta_angle)
        delta_angle = alpha_angle(mid_angle, self.endangle)
        delta2 = positive_angle(delta_angle)
        if np.isclose(delta1, delta2, rtol=rtol, atol=atol):
            result['slices'] = 1
            result['slices_half'] = 2
            result['halfslice'] = 1
            logger.debug("end of check_delta: One Area in the middle")
        else:
            result['middlelist'] = [mid_angle]
            result['slices'] = 0
            logger.debug("end of check_one_area: Area somewhere")
        return result

    def check_pairs_in_delta_list(self,
                                  delta_list,
                                  deltas,
                                  rtol=1e-3, atol=1e-2):
        if len(deltas) < 2:
            return False
        if len(delta_list) < 5:
            return False
        n1, d1 = deltas[0]
        n2, d2 = deltas[1]
        if not abs(n1 - n2) == 1:
            return False
        # check without first/last
        delta0 = delta_list[1]
        delta1 = delta_list[2]
        for delta2 in delta_list[3:-1]:
            if not np.isclose(delta0, delta2, rtol=rtol, atol=atol):
                return False
            delta0 = delta1
            delta1 = delta2
        logger.debug("** Pairs available **")
        return True

    def check_pairs_of_areas(self,
                             mid_delta_list,
                             delta,
                             geom_alpha,
                             rtol=1e-3, atol=1e-2):
        logger.debug("begin of check_pairs_of_areas")
        if len(mid_delta_list) < 2:
            return None

        logger.debug("Mid-Delta-List")
        [logger.debug(" -- mid=%s,  delta=%s",m, d) for m, d in mid_delta_list]

        # check
        mid_list = []
        m0, d0 = mid_delta_list[0]
        m1, d1 = mid_delta_list[1]
        if np.isclose(delta, d1, rtol=rtol, atol=atol):
            mx = (m0 + m1) / 2
            mid_list.append(mx)

        m0, d0 = mid_delta_list[2]
        dx = d1
        for m2, d2 in mid_delta_list[2:-1]:
            #logger.debug("compare %s and %s", d0, d2)
            if not np.isclose(d0, d2, rtol=rtol, atol=atol):
                logger.debug("end of check_pairs_of_areas: bad pairs")
                return None
            if np.isclose(delta, d2, rtol=rtol, atol=atol):
                mx = (m1 + m2) / 2
                mid_list.append(mx)
            d0 = d1
            m1 = m2
            d1 = d2

        if not mid_list:
            return None

        logger.debug("New Mids: %s", mid_list)
        delta = positive_angle(mid_list[0] * 2)
        delta_list = [delta]
        m0 = mid_list[0]
        for m1 in mid_list[1:]:
            delta = positive_angle(alpha_angle(m0, m1))
            delta_list.append(delta)
            m0 = m1

        delta_angle = alpha_angle(m1, geom_alpha)
        if np.isclose(delta_angle, np.pi*2, rtol=1e-4, atol=1e-4):
            logger.debug("Last Area is in the middle of endangle")
            delta_angle = 0.0
        delta = positive_angle(delta_angle * 2)
        delta_list.append(delta)
        logger.debug("New delta-list: %s", delta_list)
        logger.debug("end of check_pairs_of_areas")
        return delta_list

    def check_first_last_difference(self, delta_list, deltas):
        logger.debug("begin check_first_last_difference")
        if np.isclose(delta_list[0], delta_list[-1],
                      rtol=self.rtol, atol=self.atol):
            logger.debug("end check_first_last_difference: first/last equal")
            return 0
        if len(deltas) != 3:
            logger.debug("end check_first_last_difference: not 3 deltas")
            return 0
        logger.debug(">> 3 Deltas <<")
        n1, d1 = deltas[0]
        n2, d2 = deltas[1]
        n3, d3 = deltas[2]
        if not (n2 == 1 and n3 == 1):
            logger.debug("end check_first_last_difference: first/last diff")
            return 0
        dx = (d2 + d3) / 2
        if not np.isclose(dx, d1, rtol=self.rtol, atol=self.atol):
            logger.debug("end check_first_last_difference: bad deltas")
            return 0
        dx = (delta_list[0] + delta_list[-1]) / 2
        if not np.isclose(dx, d1, rtol=self.rtol, atol=self.atol):
            logger.debug("end check_first_last_difference: bad deltas")
            return 0

        logger.debug("end check_first_last_difference => %s", n1 + 1)
        return n1 + 1

    def concatenate_results(self, check_rslt):
        if len(check_rslt) < 2:
            return check_rslt
        #   ----------
        def get_slices(rslt):
            if rslt['slices'] is None:
                return 0
            if rslt.get('startdelta', 0.0) != 0.0:
                return 0
            if rslt.get('halfslice', None) is not None:
                return 0
            return rslt['slices']
        #   ------
        new_rslt = []
        size1, n1, rslt1 = check_rslt[0]
        slices1 = get_slices(rslt1)
        for size2, n2, rslt2 in check_rslt[1:]:
            slices2 = get_slices(rslt2)
            if slices1 > 0 and \
               slices1 == slices2 and \
               np.isclose(size1, size2, rtol=1e-01, atol=1e-01) and \
               (rslt1.get('delta_corr', 0.0) == 0.0 or \
                rslt2.get('delta_corr', 0.0) == 0.0):
                # concatenate
                rslt1['slices'] = None
                rslt2['delta_corr'] = 0.0
                rslt2['areas'] = rslt2['areas'] + rslt1['areas']
                logger.debug("concatenate results(size %s)", size1)
            else:
                new_rslt.append((size1, n1, rslt1))
            slices1, size1, n1, rslt1 = (slices2, size2, n2, rslt2)

        new_rslt.append((size1, n1, rslt1))
        return new_rslt

    def get_symmetry_parts(self, check_rslt):
        max_size = 0
        max_areas = 0
        max_slices = 0
        parts_possible = None
        self.delta_angle_corr = None
        unsure_sym = False

        check_rslt = [(size, n, rslt) for n, (size, rslt) in enumerate(check_rslt)
                      if rslt['slices'] is not None]
        check_rslt.sort(reverse=True)
        check_rslt = self.concatenate_results(check_rslt)
        [logger.debug("Result #%s: %s, %s", n, size, rslt) for size, n, rslt in check_rslt]

        rtol = 1e-3
        atol = 1e-2

        missing_middles = []
        halfslice = []
        start_delta = None
        start_delta_corr = 0.0

        with_angle_corr = 0
        without_angle_corr = 0
        maybe_angle_korr = 0

        for size, n, rslt in check_rslt:
            areas = rslt['areas']
            slices = rslt['slices']
            if slices is None:
                continue

            size = rslt['areasize']
            angle_corr = rslt.get('delta_corr', 0.0)

            if rslt.get('halfslice', 0) == 1:
                halfslice.append(rslt)

            if slices > 0:
                max_size = max(max_size, size)
                max_areas = max(max_areas, areas)
                max_slices = max(max_slices, slices)
                missing_middles += rslt.get('missing_middles', [])
                area = rslt['area']
                if not (np.isclose(area.min_dist,
                                   self.geom.min_radius,
                                   rtol=1e-4, atol=1e-3) and \
                        np.isclose(area.max_dist,
                                   self.geom.max_radius,
                                   rtol=1e-4, atol=1e-3)):
                    if rslt.get('delta_corr', 0.0) == 0.0:
                        without_angle_corr += areas * size
                    else:
                        with_angle_corr += areas * size
                else:
                    maybe_angle_korr += areas * size

        logger.debug("max size: %s,  max areas: %s", max_size, max_areas)
        logger.debug("Angle-Corrections: %s Yes,  %s No, %s Maybe",
                     with_angle_corr,
                     without_angle_corr,
                     maybe_angle_korr)

        if np.isclose(with_angle_corr, without_angle_corr):
            with_angle_corr = (maybe_angle_korr > 0)
        else:
            with_angle_corr = (with_angle_corr > without_angle_corr)

        #   -------------------------
        def get_halfslice_counterpart(rslt):
            if rslt.get('halfslice', 0) != 2:
                return None
            for half in halfslice:
                alpha1 = half['alpha']
                height1 = half['height']
                alpha2 = rslt['alpha'] * 2.0
                height2 = rslt['height']
                logger.debug("-- height: %s / %s", height1, height2)
                logger.debug("-- alpha: %s / %s", alpha1, alpha2)
                if np.isclose(height1, height2, rtol=rtol, atol=atol) and \
                   np.isclose(alpha1, alpha2, rtol=rtol, atol=atol):
                    return half
            return None
        #   -----------
        if halfslice:
            logger.debug("%s halfslice [1] found", len(halfslice))

            for size, n, rslt in check_rslt:
                half = get_halfslice_counterpart(rslt)
                if half:
                    logger.debug("Halfslice counterpart found")
                    alpha1 = half['alpha']
                    height1 = half['height']
                    alpha2 = rslt['alpha'] * 2.0
                    height2 = rslt['height']
                    logger.debug("-- height: %s / %s", height1, height2)
                    logger.debug("-- alpha: %s / %s", alpha1, alpha2)
                    if np.isclose(height1, height2, rtol=rtol, atol=atol) and \
                       np.isclose(alpha1, alpha2, rtol=rtol, atol=atol):
                        logger.debug("halfslice result: %s", half)
                        slices = None
                        rslt['slices'] = slices
                        half['slices'] = half['slices_half']

        angle_corr_slices = 0
        angle_corr_size = 0
        angle_corr_areas = 0
        angle_corr_airgap = False

        for size, n, rslt in check_rslt:
            areas = rslt['areas']
            slices = rslt['slices']
            size = rslt['areasize']

            if slices is None:  # ignore it
                continue

            delta_angle_corr = rslt.get('delta_corr', 0.0)
            # Angle Correction
            if with_angle_corr and self.delta_angle_corr is None:
                self.delta_angle_corr = delta_angle_corr
                angle_corr_slices = slices
                angle_corr_size = size
                angle_corr_areas = areas
                angle_corr_airgap = rslt.get('airgap', False)
            else:
                if with_angle_corr and self.delta_angle_corr != delta_angle_corr:
                    unsure_sym = True
                    if slices > angle_corr_slices and \
                       size > angle_corr_size * 0.5:
                        logger.debug("Angle Correction")
                        self.delta_angle_corr = delta_angle_corr
                        angle_corr_slices = slices
                        angle_corr_size = size
                        angle_corr_areas = areas
                        angle_corr_airgap = rslt.get('airgap', False)
                    elif angle_corr_airgap and \
                         slices >= angle_corr_slices and \
                         areas > angle_corr_areas and \
                         size > angle_corr_size * 0.2:
                        logger.debug("Angle Correction")
                        self.delta_angle_corr = delta_angle_corr
                        angle_corr_slices = slices
                        angle_corr_size = size
                        angle_corr_areas = areas
                        angle_corr_airgap = rslt.get('airgap', False)
            if slices:
                if start_delta is None:
                    start_delta = rslt['startdelta']
                    start_delta_corr = start_delta
                elif not rslt['startdelta'] != 0.0:
                    if start_delta_corr == 0.0:
                        start_delta_corr = rslt['startdelta']
                    elif not np.isclose(rslt['startdelta'], start_delta_corr,
                                        rtol=rtol, atol=atol):
                        slices = 0  # bad
                        rslt['slices'] = slices

            if slices is not None and slices <= 1:  # critical
                if areas <= max(1, max_areas / 5):
                    if size < max_size / 25:
                        slices = None
            if slices == 0:
                middles = rslt.get("middlelist", [])
                if self.check_missing_areas(middles,
                                            missing_middles):
                    logger.debug("Symmetry-Destroyer destroyed")
                    slices = None

            if slices == 1:
                # symmetry killer
                if areas < max(2, max_areas / 6):
                    if size < max_size * 0.05:
                        slices = None  # ignore tiny areas


            parts_possible = self.calc_parts(parts_possible, slices)

        if unsure_sym:
            logger.warning("Warning: unsure symmetry")

        if parts_possible is None:
            parts_possible = 0
        return parts_possible, start_delta

    def check_missing_areas(self, middles, missing_middles):
        logger.debug("check_missing_areas")
        logger.debug(" -- mids = %s", middles)
        logger.debug(" -- missing mids = %s", missing_middles)
        if not missing_middles:
            return False
        if len(middles) == 0 or len(middles) > 2:
            return False

        for m in middles:
            mlist = [mm for mm, i in missing_middles
                     if np.isclose(m, mm, rtol=1e-3, atol=1e-3)]
            if not mlist:
                return False
        return True

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

    def create_cut_lines(self, parts, start_delta):
        self.geom.clear_cut_lines()
        for alpha in self.symmetry_lines(parts,
                                         self.startangle,
                                         start_delta,
                                         self.endangle):
            plus = self.geom.max_radius / 10
            min_radius = max(10, self.geom.min_radius - plus)
            p1 = point(self.geom.center, min_radius, alpha)
            p2 = point(self.geom.center, self.geom.max_radius + plus, alpha)
            line = Line(Element(start=p1, end=p2))
            line.init_attributes(color='green')
            self.geom.add_cut_line(line)

    def symmetry_lines(self, parts, startangle, start_delta, endangle):
        logger.debug("begin symmetry_lines from %s to %s with start %s",
                     startangle,
                     endangle,
                     start_delta)

        if less_equal(endangle, startangle):
            endangle += 2*np.pi

        sym = self.geom_part * parts
        delta = 2*np.pi/sym
        if start_delta == 0.0:
            if not is_same_angle(startangle, endangle):
                start_delta = delta

        sym_startangle = startangle + start_delta
        start = startangle + start_delta
        while less(start, endangle):
            yield start
            start += delta

        if is_same_angle(startangle, endangle):
            yield start

        # Damit man anschliessend ohne Umstände schneiden kann.
        self.set_symmetry_parameters(sym_startangle, parts, delta)
        logger.debug("end symmetry_lines")

    def set_symmetry_parameters(self, startangle, parts, delta):
        self.geom.sym_startangle = startangle
        self.geom.sym_endangle = startangle + delta
        self.geom.sym_slices = parts
        self.geom.sym_slice_angle = delta
        self.geom.sym_area = Area([], (0,0), 0.0)
        self.geom.sym_area.sym_startangle = self.geom.sym_startangle
        self.geom.sym_area.sym_endangle = self.geom.sym_endangle

    def check_symmetry_of_mirror(self, mirror_geom, mirrorangle):
        logger.debug("begin of Symmetry::check_symmetry_of_mirror")
        assert(mirror_geom is not None)

        axis_p = point(self.geom.center, self.geom.max_radius, mirrorangle)
        axis_m = line_m(self.geom.center, axis_p)
        axis_n = line_n(self.geom.center, axis_m)

        def counterpart_found(node, mirror_node, nodes, rtol, atol):
            hit_sloppy = 0
            for n in nodes:
                if points_are_close(mirror_node, n, rtol, atol):
                    logger.debug(" ---- %s is %s", node, mirror_node)
                    return 1, 1
                if points_are_close(round_point(mirror_node, 1),
                                    round_point(n, 1),
                                    rtol, atol):
                    logger.debug(" ++++ %s is %s", node, mirror_node)
                    hit_sloppy = 1

            logger.debug(" >>>> %s is NOT %s  (%s)",
                         node, mirror_node, hit_sloppy)
            return 0, hit_sloppy

        def check_differences(geom, mirror_geom):
            geom_ag_nodes = []
            geom_nodes = [n for n in geom.g.nodes() if not (n in geom_ag_nodes)]

            hits = 0
            hits_sloppy = 0
            for n in geom_nodes:
                mirror_n = mirror_point(n, geom.center, axis_m, axis_n)
                hit, hit_sloppy = counterpart_found(n, mirror_n,
                                                    mirror_geom.g.nodes(),
                                                    # self.rtol,
                                                    1e-3,
                                                    # self.atol):
                                                    1e-2)
                hits += hit
                hits_sloppy += hit_sloppy

            min_nodes = min(len(geom_nodes), int(len(geom_nodes) * 0.95) + 1)
            logger.debug("Nodes=%s,  Counterparts=%s (sloppy=%s)",
                         len(geom_nodes), hits, hits_sloppy)
            if hits < min_nodes:
                f = hits / len(geom_nodes)
            else:
                f = 1.0
            if hits_sloppy < min_nodes:
                f_sloppy = hits_sloppy / len(geom_nodes)
            else:
                f_sloppy = 1.0
            return f, f_sloppy

        # ----------------
        logger.debug("check geom - mirror")
        f1, f1_sloppy = check_differences(self.geom, mirror_geom)
        logger.debug("check mirror - geom")
        f2, f2_sloppy = check_differences(mirror_geom, self.geom)
        logger.debug("Factor 1: %s,  2: %s", f1, f2)
        if f1 >= 0.99 or f2 >= 0.99:
            ok = not (f1 < 0.9 or f2 < 0.9)
        else:
            ok = f1 > 0.97 and f2 > 0.97
        if not ok:
            if f1_sloppy > 0.97 and f2_sloppy > 0.97:
                logger.debug("  (A sloppy mirror found, but ignored)")
                ok = True
        logger.debug("end of Symmetry::check_symmetry_of_mirror => %s", ok)
        return ok
