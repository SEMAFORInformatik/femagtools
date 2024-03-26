# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.concat
  ~~~~~~~~~~~~~~~~~~~~~~~

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function

import numpy as np
import logging
import sys
from femagtools.dxfsl.shape import Element, Shape, Circle, Arc, Line
from femagtools.dxfsl.shape import is_Circle, is_Arc, is_Line
from femagtools.dxfsl.functions import points_are_close, positive_angle, elevation_angle
from femagtools.dxfsl.functions import Timer, alpha_line
from femagtools.dxfsl.journal import getJournal
import io
import time

logger = logging.getLogger('femagtools.concat')


#############################
#           concat          #
#############################


class Concatenation(object):
    def __init__(self,
                 rtol=1e-04,
                 atol=1e-04,
                 ndec=6,
                 mdec=1):
        self.rtol = rtol
        self.atol = atol
        self.ndec = ndec
        self.mdec = mdec
        self.line_y0 = Line(Element(start=(-1e+8, 0.0), end=(1e+8, 0.0)))
        self.line_x0 = Line(Element(start=(0.0, -1e+8), end=(0.0, 1e+8)))
        self.journal = getJournal()

        logger.debug("Concatenation(rtol=%s, atol=%s, ndec=%s)", rtol, atol, ndec)

    def __str__(self):
        return "rtol: {}\n".format(self.rtol) + \
               "atol: {}\n".format(self.atol) + \
               "ndec: {}\n".format(self.ndec)

    def concatenate_line_elements(self,
                                  src_elements,
                                  dest_elements):
        if len(src_elements) < 2:
            return 0

        src_list = []
        for e in src_elements:
            if not e.has_attribute("del"):
                e_p1, ep_2 = e.points_sorted()
                src_list.append((e_p1, e))
        if len(src_list) < 2:
            return 0

        logger.debug("begin of concatenate_line_elements (%s elements)",
                     len(src_list))

        src_list.sort()
        e1_p1, e1 = src_list[0]
        new_elements = []
        count = 0
        for e2_p1, e2 in src_list[1:]:
            el_new = e1.concatenate(None, None, e2,
                                    rtol=1e-06,
                                    atol=1e-06,
                                    mdec=self.mdec,
                                    overlapping=True)
            if el_new:
                logger.debug("NEW %s from %s to %s",
                             el_new.classname(),
                             el_new.p1,
                             el_new.p2)
                el_new.set_attribute("new")
                e1.set_attribute("del")
                e2.set_attribute("del")
                logger.debug("OLD %s from %s to %s",
                             e1.classname(),
                             e1.p1,
                             e1.p2)
                logger.debug("OLD %s from %s to %s",
                             e2.classname(),
                             e2.p1,
                             e2.p2)
                el_p1, el_p2 = el_new.points_sorted()
                dest_elements.append(el_new)
                src_elements.append(el_new)
                e1 = el_new
                count += 1
            else:
                e1 = e2

        logger.debug("end of concatenate_line_elements: %s concatenations", count)
        return count

    def concatenate_arc_elements(self,
                                 src_elements,
                                 dest_elements):
        if len(src_elements) < 2:
            return 0
        
        src_list = [e for e in src_elements if not e.has_attribute("del")]
        if len(src_list) < 2:
            return 0

        e1 = src_list[0]
        new_elements = []
        count = 0
        for e2 in src_list[1:]:
            el_new = e1.concatenate(None, None, e2,
                                    rtol=1e-06,
                                    atol=1e-06,
                                    mdec=self.mdec,
                                    overlapping=True)
            if el_new:
                logger.debug("NEW %s from %s to %s",
                             el_new.classname(),
                             el_new.p1,
                             el_new.p2)
                el_new.set_attribute("new")
                e1.set_attribute("del")
                e2.set_attribute("del")
                logger.debug("OLD %s from %s to %s",
                             e1.classname(),
                             e1.p1,
                             e1.p2)
                logger.debug("OLD %s from %s to %s",
                             e2.classname(),
                             e2.p1,
                             e2.p2)
                dest_elements.append(el_new)
                src_elements.append(el_new)
                e1 = el_new
                count += 1
            else:
                e1 = e2
        return count

    def concatenate_matching_line_elements(self,
                                           elements,
                                           only_tiny=False):
        logger.debug("Begin of concatenate_matching_line_elements")
        timer = Timer(start_it=True)

        if only_tiny:
            tiny = [e for e in elements if is_Line(e) and e.has_attribute("tiny")]
            if not tiny:
                timer.stop("-- no tiny lines found in %0.4f seconds --")
                return 0  # no tiny line exists
        
        def axis_point(e, dec=1):
            if np.isclose(e.p1[0], e.p2[0]):
                return (round(e.p1[0], dec), 0.0)
            if np.isclose(e.p1[1], e.p2[1]):
                return (0.0, round(e.p1[1], dec))
            px = e.intersect_line(self.line_y0, all=True)
            py = e.intersect_line(self.line_x0, all=True)
            if not (px or py):
                logger.error("axis_point: FATAL ERROR !!")
                return (0.0, 0.0)

            if not px:
                py = round(py[0][1], dec)
                return (0.0, py)
            if not py:
                px = round(px[0][0], dec)
                return (px, 0.0)
            
            px = round(px[0][0], dec)
            py = round(py[0][1], dec)
            
            if abs(px) > abs(py):
                return (0.0, py)
            return (px, 0.0)

        def elevation(e):
            return round(elevation_angle(alpha_line(e.p1, e.p2)), 2)

        def axis_points_are_close(a1, a2):
            #return a1 == a2
            if a1[0]:
                d = abs(a1[0] - a2[0])  # x
            else:
                d = abs(a1[1] - a2[1])  # y
            if d < 0.2:
                logger.debug("axis_points_are_close: %s - %s = %s", a1, a2, d)
            return d < 0.2

        elmts = [(round(e.m(999999.0), self.mdec),
                  axis_point(e, dec=4),
                  e.p1,
                  e) for e in elements if is_Line(e)]
        if not elmts:
            return 0  # no lines available

        elmts.sort()

        logger.debug("Concatenate Line Elements")
        #for m, a, p, e in elmts:
        #    logger.debug("Line from %s to %s [m=%s, a=%s]",
        #                 e.p1, e.p2, m, a)
        logger.debug("*************************")

        lines_available = len(elmts)
        count = 0
        if len(elmts) < 2:
            timer.stop()
            return 0

        new_elements = []
        m1, a1, p1, e1 = elmts[0]
        line_elements = [e1]
        for m2, a2, p2, e2 in elmts[1:]:
            if m1 == m2 and axis_points_are_close(a1, a2):
                line_elements.append(e2)
            else:
                if len(line_elements) > 1:
                    logger.debug("Begin Concatenate %s Lines: m=%s, a=%s",
                                 len(line_elements), m1, a1)
                    for l in line_elements:
                        logger.debug("-- %s", l)
                    line_count = self.concatenate_line_elements(line_elements, elements)
                    while line_count > 0:
                        count += line_count
                        line_count = self.concatenate_line_elements(line_elements, elements)
                    logger.debug("End Concatenate Lines: m=%s, a=%s", m1, a1)

                m1 = m2
                a1 = a2
                line_elements = [e2]

        logger.debug("Begin Concatenate Lines: m=%s, a=%s", m1, a1)
        line_count = self.concatenate_line_elements(line_elements, elements)
        while line_count > 0:
            count += line_count
            line_count = self.concatenate_line_elements(line_elements, elements)
        logger.debug("End Concatenate Lines: m=%s, a=%s", m1, a1)

        timer.stop("-- {} lines concatenated in %0.4f seconds --".format(count))
        line_list = [e for e in elements
                    if not e.has_attribute("del") and is_Line(e)]
        logger.debug("== %s of %s Line-Elements are remaining ==",
                     len(line_list),
                     lines_available)
        logger.debug("End of concatenate_matching_line_elements")
        return count

    def concatenate_matching_arc_elements(self,
                                          elements,
                                          only_tiny=False):
        logger.debug("Begin of concatenate_matching_arc_elements")
        timer = Timer(start_it=True)

        if only_tiny:
            tiny = [e for e in elements if is_Arc(e) and e.has_attribute("tiny")]
            if not tiny:
                timer.stop("-- no tiny arcs found in %0.4f seconds --")
                logger.debug("End of concatenate_matching_arc_elements")
                return 0  # no tiny arc exists

        def center_rounded(c):
            return (round(c[0], 6),round(c[1], 6))

        elmts = [(center_rounded(e.center),
                  round(e.radius, 6),
                  positive_angle(e.startangle),
                  e) for e in elements if is_Arc(e) or is_Circle(e)]
        if not elmts:
            timer.stop("-- no arcs found in %0.4f seconds --")
            logger.debug("End of concatenate_matching_arc_elements")
            return 0  # no arcs available

        elmts.sort()

        logger.debug("Concatenate Arc Elements")
        #for c, r, a, e in elmts:
        #    logger.debug("%s from %s to %s [c=%s, r=%s, a=%s]",
        #                 e.classname(), e.p1, e.p2, c, r, a)
        logger.debug("*************************")
 
        arcs_available = len(elmts)
        count = 0

        c1, r1, a1, e1 = elmts[0]
        arc_elements = [e1]
        for c2, r2, a2, e2 in elmts[1:]:
            if points_are_close(c1, c2) and np.isclose(r1, r2):
                arc_elements.append(e2)
            else:
                logger.debug("Begin Concatenate %s Arcs: c=%s, r=%s",
                             len(arc_elements), c1, r1)
                arc_count = self.concatenate_arc_elements(arc_elements, elements)
                while arc_count > 0:
                    count += arc_count
                    arc_count = self.concatenate_arc_elements(arc_elements, elements)
                logger.debug("End Concatenate Arcs: c=%s, r=%s", c1, r1)

                c1 = c2
                r1 = r2
                arc_elements = [e2]

        logger.debug("Begin Concatenate Arcs: c=%s, r=%s", c1, r1)
        arc_count = self.concatenate_arc_elements(arc_elements, elements)
        while arc_count > 0:
            count += arc_count
            arc_count = self.concatenate_arc_elements(arc_elements, elements)
        logger.debug("End Concatenate Arcs: c=%s, r=%s", c1, r1)

        timer.stop("-- {} arcs concatenated in %0.4f seconds --".format(count))

        arc_list = [e for e in elements
                    if not e.has_attribute("del") and
                    (is_Arc(e) or is_Circle(e))]
        logger.debug("== %s of %s Arc-Elements are remaining ==",
                     len(arc_list),
                     arcs_available)
        logger.debug("End of concatenate_matching_arc_elements")
        return count

    def concatenate_matching_elements(self,
                                      src_elements,
                                      main=False,
                                      only_tiny=False):
        logger.debug("Begin of concatenate_matching_elements")
        timer = Timer(start_it=True)
        count_lines = count_arcs = 0
        count_lines = self.concatenate_matching_line_elements(src_elements,
                                                              only_tiny=only_tiny)
        count_arcs = self.concatenate_matching_arc_elements(src_elements,
                                                            only_tiny=only_tiny)
        if only_tiny:
            if count_lines == 0 and count_arcs == 0:
                #  no action
                timer.stop("-- no concatenations in %0.4f seconds --")
                return False, src_elements

        count = count_lines + count_arcs
        if main:
            self.journal.put_concat_lines(count_lines)
            self.journal.put_concat_arcs(count_arcs)

        new_list = [e for e in src_elements if not e.has_attribute("del")]

        t = timer.stop("-- {} concatenations in %0.4f seconds --".format(count))
        if main:
            self.journal.put('time_concatenation', t)

        logger.debug("End of concatenate_matching_elements")
        return count, new_list
