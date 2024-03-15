# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.areabuilder
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function

import numpy as np
import logging
import sys
from femagtools.dxfsl.shape import Element, Shape, Circle, Arc, Line
from femagtools.dxfsl.shape import is_Circle, is_Arc, is_Line
from femagtools.dxfsl.area import Area
from femagtools.dxfsl.functions import points_are_close, nodes_are_equal, distance
from femagtools.dxfsl.functions import normalise_angle, positive_angle, point
from femagtools.dxfsl.functions import alpha_line, alpha_points, alpha_angle
from femagtools.dxfsl.functions import less, is_same_angle
from femagtools.dxfsl.functions import Timer
from femagtools.dxfsl.journal import getJournal
import io
import time

logger = logging.getLogger('femagtools.concat')


def log_lefthand(left, edge1, edge2):
    if left:
        edge1.log_edge("**** left")
    else:
        edge2.log_edge("**** left")


#############################
#        areabuilder        #
#############################


class EdgeInfo(object):
    def __init__(self,
                 edge_data=None,
                 n1 = None,
                 n2 = None,
                 handmade=False):
        assert(edge_data)
        assert(n1)
        assert(n2)

        if not handmade:
            e = edge_data.get('object', None)
            if not e:
                raise ValueError("Fatal: no object in edge_data %s", edge_data)          
            x = e.get_node_number(n1)
        else:
            e = None
            x = 0
            
        self.n1 = n1
        self.n2 = n2
        self.data = edge_data
        self.element = e
        self.x = x
        self.tracked = edge_data.get(x, False)
        self.alpha = 0.0
        self.errors = 0
        self.startangle = None
        self.angle = None
        self.name = None

    def log_edge(self, text):
        logger.debug("%s: %s[%s] from %s to %s [%s]",
                     text,
                     self.classname(),
                     self.name,
                     self.n1, self.n2,
                     self.angle)

    def log_name(self, name):
        self.name = name

    def get_reverse_edge(self):
        return EdgeInfo(edge_data=self.data, n1=self.n2, n2=self.n1)

    def is_circle(self):
        return is_Circle(self.element)

    def is_arc(self):
        return is_Arc(self.element)

    def is_line(self):
        return is_Line(self.element)

    def classname(self):
        return self.element.classname()

    def is_start_edge(self):
        return self.startangle is not None

    def set_tracked(self):
        self.data[self.x] = True  # footprint
        self.tracked = True

    def set_all_tracked(self):
        self.data[1] = True  # footprint
        self.data[2] = True  # footprint
        self.tracked = True

    def is_tracked(self):
        return self.tracked

    def set_start_angle(self):
        start_angle = positive_angle(self.n2_angle_ingoing())
        self.startangle = start_angle

    def set_direction_angle(self, start_angle):
        angle = positive_angle(alpha_angle(start_angle, self.n1_angle_ingoing()))
        if is_same_angle(0.0, angle, atol=0.008):  # 1/2 degree
            # reverse direction
            if self.n1_direction_righthand():
                self.angle = 2.0 * np.pi
            elif self.n1_direction_lefthand():
                self.angle = 0.0
            else:
                self.angle = angle
        else:
            self.angle = angle

    def n1_angle_outgoing(self):
        return self.element.get_alpha(self.n1)

    def n1_angle_ingoing(self):
        return normalise_angle(self.element.get_alpha(self.n1) + np.pi)

    def n2_angle_outgoing(self):
        return self.element.get_alpha(self.n2)

    def n2_angle_ingoing(self):
        return normalise_angle(self.element.get_alpha(self.n2) + np.pi)

    def n1_direction_lefthand(self):
        if not self.is_arc():
            return False
        alpha_in = self.n1_angle_ingoing()
        alpha_c = alpha_line(self.n1, self.element.center)
        angle = positive_angle(alpha_angle(alpha_in, alpha_c))
        return angle < np.pi

    def n1_direction_righthand(self):
        return not self.n1_direction_lefthand()

    def radius(self):
        if self.is_arc():
            return self.element.radius
        return 0.0

    def line_arc_close_together(self, info):
        assert(self.is_arc())
        assert(info.is_line())

        #if self.is_start_edge() or info.is_start_edge():
        #    return False

        dist_n2 = distance(self.element.center, info.n2)
        return np.isclose(self.element.radius, dist_n2, 0.005, 0.005)

    def myself_direction_lefthand(self, start_edge, nbr_edge, builder, ignore_start=False):
        logger.debug("start of myself_direction_lefthand")
        self.log_edge("---> self")
        nbr_edge.log_edge("---> nbr")

        myself_angle = self.angle
        other_angle = nbr_edge.angle

        if is_same_angle(0.0, myself_angle):
            # 360 or 0 degrees => turn 180 degrees
            logger.debug("-- ATTENTION: myself %s turns nearly 180 degrees", self.classname())
            logger.debug("   the angle is %s", myself_angle)
            if start_edge.is_arc():
                if self.is_line() and not ignore_start:
                    if start_edge.line_arc_close_together(self):
                        logger.debug("START ARC and SELF LINE close")
                        arc_edge = start_edge.get_reverse_edge()  # reverse start edge
                        arc_edge.log_name("REVERSE")
                        arc_edge.set_direction_angle(start_edge.startangle)
                        left = not arc_edge.arc_line_direction_lefthand(start_edge, self, builder)
                        log_lefthand(left, self, arc_edge)
                        if left:  # self is left
                            myself_angle = 0.0
                        else:
                            myself_angle = np.pi * 2.0
                        left = myself_angle > other_angle

                        log_lefthand(left, self, nbr_edge)
                        logger.debug("#3: end of myself_direction_lefthand: ==> %s", left)
                        return left
                    
                if start_edge.n1_direction_lefthand():
                    myself_angle = 0.0
                else:
                    myself_angle = np.pi * 2.0
    
            elif self.is_arc():
                if not ignore_start:
                    if self.line_arc_close_together(start_edge):
                        logger.debug("START LINE and SELF ARC close")
                        line_edge = start_edge.get_reverse_edge()  # reverse start edge
                        line_edge.log_name("REVERSE")
                        line_edge.set_direction_angle(start_edge.startangle)
                        left = self.arc_line_direction_lefthand(start_edge, line_edge, builder)
                        log_lefthand(left, self, line_edge)
                        if left:  # self is left
                            myself_angle = 0.0
                        else:
                            myself_angle = np.pi * 2.0
                        left = myself_angle > other_angle

                        log_lefthand(left, self, nbr_edge)
                        logger.debug("#4: end of myself_direction_lefthand: ==> %s", left)
                        return left

                if self.n1_direction_lefthand():
                    myself_angle = 0.0
                else:
                    myself_angle = np.pi * 2.0

        if is_same_angle(0.0, other_angle, atol=0.008):  # 1/2 degree
            # 360 or 0 degrees => turn 180 degrees
            logger.debug("-- ATTENTION: other %s turns nearly 180 degrees", nbr_edge.classname())
            logger.debug("   the angle is %s", other_angle)
            if start_edge.is_arc():
                if nbr_edge.is_line() and not ignore_start:
                    if start_edge.line_arc_close_together(nbr_edge):
                        logger.debug("START ARC and SELF LINE close")
                        arc_edge = start_edge.get_reverse_edge()  # reverse start edge
                        arc_edge.set_direction_angle(start_edge.startangle)
                        left = arc_edge.arc_line_direction_lefthand(start_edge, nbr_edge, builder)
                        log_lefthand(left, self, line_edge)
                        logger.debug("#5: end of myself_direction_lefthand: ==> %s", not left)
                        return not left

                if start_edge.n1_direction_lefthand():
                    other_angle = 0.0
                else:
                    other_angle = np.pi * 2.0

            elif nbr_edge.is_arc():
                if not ignore_start:
                    if nbr_edge.line_arc_close_together(start_edge):
                        logger.debug("START LINE and NEIGHBOR ARC close")
                        line_edge = start_edge.get_reverse_edge()  # reverse start edge
                        line_edge.set_direction_angle(start_edge.startangle)
                        left = nbr_edge.arc_line_direction_lefthand(start_edge, line_edge, builder)
                        logger.debug("#6: end of myself_direction_lefthand: ==> %s", left)
                        return left

                if nbr_edge.n1_direction_lefthand():
                    other_angle = 0.0
                else:
                    other_angle = np.pi * 2.0

        logger.debug("-- angles: myself = %s,  other = %s",
                     myself_angle, other_angle)
        if not is_same_angle(myself_angle, other_angle, atol=0.008):  # 1/2 degree
            logger.debug("-- angles are different")
            left = myself_angle > other_angle
            log_lefthand(left, self, nbr_edge)
            logger.debug("end of myself_direction_lefthand: ==> %s", left)
            return left

        logger.debug("nearly same direction")

        if self.is_line():
            if nbr_edge.is_line():
                return self.angle > nbr_edge.angle
            left = nbr_edge.arc_line_direction_lefthand(start_edge, self, builder)
            log_lefthand(left, nbr_edge, self)
            logger.debug("#1: end of myself_direction_lefthand: ==> %s", not left)
            return not left

        if nbr_edge.is_line():
            left = self.arc_line_direction_lefthand(start_edge, nbr_edge, builder)
            log_lefthand(left, self, nbr_edge)
            logger.debug("#2: end of myself_direction_lefthand: ==> %s", left)
            return left

        # two arcs
        left = self.myself_lefthand_side(nbr_edge)
        log_lefthand(left, self, nbr_edge)
        logger.debug("end of myself_direction_lefthand: ==> %s", left)
        return left

    def myself_lefthand_side(self, edge):
        # two edges with same direction
        if self.is_arc():
            if edge.is_line():
                return self.n1_direction_lefthand()
            
            edge_left = edge.n1_direction_lefthand()
            if self.n1_direction_lefthand():
                if not edge_left:
                    return True
                return self.radius() < edge.radius()
            else:
                if edge_left:
                    return False
                return self.radius() > edge.radius()
        else:
            if edge.is_arc():
                return edge.n1_direction_righthand()
        return True  # two lines

    def arc_line_direction_lefthand(self, start_edge, line_edge, builder):
        logger.debug("begin of arc_line_direction_lefthand")
        assert(self.is_arc())
        assert(line_edge.is_line())
        start_edge.log_edge(":::: start")
        self.log_edge(".... self")
        line_edge.log_edge(".... line")

        if not self.line_arc_close_together(line_edge):
            logger.debug("-- line and arc not close together")
            left = self.n1_direction_lefthand()  # ARC
            logger.debug("end of arc_line_direction_lefthand: arc lefthand = %s", left)
            return left

        logger.debug("-- ATTENTION: line look ahead")
        line_edge.startangle = start_edge.startangle
        # Line is start-edge now
        next_edge = builder.next_edge_lefthand_side(line_edge)  # next of line
        next_edge.set_direction_angle(start_edge.startangle)
        left = self.myself_direction_lefthand(start_edge,
                                              next_edge,
                                              builder,
                                              ignore_start=True)
        logger.debug("end of arc_line_direction_lefthand: arc lefthand = %s", left)
        return left


class AreaBuilder(object):
    def __init__(self,
                 geom=None,
                 rtol=1e-04,
                 atol=1e-04):
        self.rtol = rtol
        self.atol = atol
        self.geom = geom
        self.area_list = []
        self.journal = getJournal()

    def __str__(self):
        return "rtol: {}\n".format(self.rtol) + \
               "atol: {}\n".format(self.atol)

    def all_neighbors(self, n):
        for n in self.geom.g.neighbors(n):
            yield n

    def set_edge_attributes(self):
        self.geom.set_edge_attributes()
    
    def create_list_of_areas(self, main=False):
        logger.debug("Begin of create_list_of_areas")
        assert(len(self.area_list) == 0)
        self.errors = 0

        def append(area_list, a):
            for area in area_list:
                if area.is_identical(a):
                    return
            area_list.append(a)

        timer = Timer(start_it=True)

        self.set_edge_attributes()
        
        for n in self.geom.nodes():
            finished = False
            while not finished:
                finished = True
                nbrs = [nbr for nbr in self.geom.g.neighbors(n)]
                for next_n in nbrs:
                    result = self.get_new_area(n, next_n)
                    if result['ok']:
                        area = result['area']
                        a = Area(area, self.geom.center, 0.0)
                        logger.debug("Area %s found", a.identifier())
                        append(self.area_list, a)

        t = timer.stop("{} areas created in %0.4f seconds".format(len(self.area_list)))
        if main:
            self.journal.put('time_build_areas', t)
        if self.errors > 0:
            logger.warning("WARNING: %s errors while creating %s areas",
                           self.errors, len(self.area_list))
        self.journal.put_areas(len(self.area_list))
        self.journal.put_area_errors(self.errors)
        logger.debug("End of create_list_of_areas: %s areas created",
                     len(self.area_list))

    def get_neighbors(self, info):
        nbrs = [n for n in self.all_neighbors(info.n2)
                if not nodes_are_equal(n, info.n1)]
        if len(nbrs) == 0:
            logger.error("FATAL ERROR: no neighbors for node %s available ???",
                         info.n2)
            return []

        neighbors = []
        for n in nbrs:
            neighbors.append(self.get_edge_info(info.n2, n))
        return neighbors
    
    def get_edge_info(self, n1, n2):
        edge_data = self.geom.g.get_edge_data(n1, n2)
        if not edge_data:
            raise ValueError("Fatal: no edge-data found from {} to {}"
                             .format(n1, n2))

        return EdgeInfo(edge_data, n1, n2)

    def get_new_area(self, start_n1, start_n2):
        info_curr = self.get_edge_info(start_n1, start_n2)
        result = {'area': None,
                  'elements': 1,
                  'msg': "<undefined>",
                  'reverse': False,
                  'ok': False}
 
        if info_curr.is_tracked():
            result['msg'] = ("<== area already tracked (%s) ***",
                             info_curr.x)
            return result

        area = [info_curr.element]
        result['area'] = area

        logger.debug('begin of get_new_area')

        if info_curr.is_circle():
            result['msg'] = "area is a circle !!"
            logger.debug("<== %s", result['msg'])
            info_curr.set_all_tracked()
            result['ok'] = True
            return result

        logger.debug("***** EDGE %s *****", 1)
        info_curr.set_start_angle()
        info_curr.log_name("START")
        info_next = self.next_edge_lefthand_side(info_curr)
        if not info_next:
            result['msg'] = ("dead end ({}, {})"
                             .format(info_curr.n1, info_curr.n2))
            logger.debug("<== %s", result['msg'])
            self.errors += 1
            return result

        prev_n1 = info_curr.n1
        next_n1 = info_next.n1
        next_n2 = info_next.n2

        alpha = normalise_angle(alpha_points(prev_n1,
                                             next_n1,
                                             next_n2))

        c = 1
        while not (nodes_are_equal(next_n1, start_n1) and
                   nodes_are_equal(next_n2, start_n2)):
            c += 1
            if c > self.geom.num_edges * 2:
                logger.error("FATAL: *** over %s elements in area ? ***",
                             self.geom.num_edges)
                sys.exit(1)

            area.append(info_next.element)

            if info_next.is_tracked():
                result['msg'] = ("FATAL: area already tracked ({}) ***"
                                 .format(info_next.x))
                logger.debug("<== %s", result['msg']) # ERROR
                logger.debug("++ %s from %s to %s",
                             info_next.classname(), info_next.n1, info_next.n2)
                result['elements'] = c
                self.errors += 1
                return result

            info_next.set_tracked()

            info_curr = info_next
            logger.debug("***** EDGE %s *****", c)
            info_next.set_start_angle()
            info_curr.log_name("START")
            info_next = self.next_edge_lefthand_side(info_curr)
            if not info_next:
                result['msg'] = ("<== dead end ({},{})"
                                 .format(info_curr.n1, info_curr.n2))
                logger.debug("<== %s", result['msg'])
                result['elements'] = c
                self.errors += 1
                return result

            prev_n1 = info_curr.n1
            next_n1 = info_next.n1
            next_n2 = info_next.n2

            a = normalise_angle(alpha_points(prev_n1,
                                             next_n1,
                                             next_n2))
            alpha += a

        logger.debug("  END OF get_new_area")

        if less(alpha, 0.0):
            result['msg'] = ("turn left expected, but it turned right ({})"
                             .format(alpha))
            logger.debug("<== %s", result['msg'])
            result['elements'] = c
            result['reverse'] = True
            # it's not an error
            return result

        result['msg'] = "area found !!"
        logger.debug("<== %s", result['msg'])
        result['elements'] = c
        result['ok'] = True
        return result

    def next_edge_lefthand_side(self, start_edge):  # current
        logger.debug("begin of next_edge_lefthand_side")
        logger.debug("-- start is %s from %s to %s",
                     start_edge.classname(),
                     start_edge.n1,
                     start_edge.n2)

        nbrs = self.get_neighbors(start_edge)
        if len(nbrs) == 0:
            logger.debug("WARNING: no neighbors available ???")
            return None  # unexpected end => appendix

        logger.debug("-- NODE %s has %s neighbors", start_edge.n2, len(nbrs))
        logger.debug("-- startangle is %s", start_edge.startangle)

        if len(nbrs) == 1:
            logger.debug("end of next_edge_lefthand_side: one neighbor")
            return nbrs[0]

        all_nbrs = []
        logger.debug("candidates are:") 
        for i, nbr in enumerate(nbrs):
            nbr.set_direction_angle(start_edge.startangle)
            nbr.log_name("EDG-{}".format(i))
            nbr.log_edge("++++ candidate")
            all_nbrs.append((nbr.angle, i, nbr))
        all_nbrs.sort()

        a1, i1, nbr1 = all_nbrs[0]  # lefthandside
        nbr1.log_edge("==== first")
        for a2, i2, nbr2 in all_nbrs[1:]:
            nbr2.log_edge("==== next")
            if nbr2.myself_direction_lefthand(start_edge, nbr1, self):
                nbr1 = nbr2
                logger.debug("-- is new lefthand neighbor")

        nbr1.log_edge(">>>> LEFT")
        logger.debug("end of next_edge_lefthand_side")
        return nbr1
