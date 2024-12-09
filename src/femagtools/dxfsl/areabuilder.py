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
from femagtools.dxfsl.area import Area, TYPE_AIR
from femagtools.dxfsl.functions import points_are_close, nodes_are_equal, distance
from femagtools.dxfsl.functions import normalise_angle, positive_angle, point
from femagtools.dxfsl.functions import alpha_line, alpha_points, alpha_angle
from femagtools.dxfsl.functions import less, is_same_angle
from femagtools.dxfsl.functions import Timer
from femagtools.dxfsl.journal import getJournal
import io
import time

logger = logging.getLogger('femagtools.areabuilder')
original_log_level = None

def disable_logging():
    logger
    global original_log_level
    # Store the current log level to restore it later
    original_log_level = logger.getEffectiveLevel()
    logger.debug("Logging level %s disabled", original_log_level)

    # Set the log level to a higher level, e.g., WARNING or CRITICAL
    logging.disable(logging.ERROR)


def enable_logging():
    #if original_log_level is not None:
    # Restore the original log level after the tests
    logging.disable(logging.NOTSET)
    logger.debug("Logging level %s enabled", original_log_level)

    
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

    def set_direction_angle(self, start_edge):
        start_angle = start_edge.startangle
        logger.debug("begin set_direction_angle: of %s", self.name)
        self.angle = positive_angle(alpha_angle(start_angle, self.n1_angle_ingoing()))
        logger.debug("set_direction_angle: angle is %s", self.angle)

        if is_same_angle(0.0, self.angle, atol=0.008):  # 1/2 degree
            # reverse direction
            logger.debug("set_direction_angle: reverse direction( nearly 180 degrees)")

            if start_edge.is_arc():
                start_lefthand = start_edge.n1_direction_lefthand()
                if self.is_arc():
                    myself_lefthand = self.n1_direction_lefthand()
                    if start_lefthand:
                        if not myself_lefthand:  # righthand => same side
                            left = self.radius() < start_edge.radius()
                        else:
                            left = False
                    else:
                        if myself_lefthand:  # lefthand => same side
                            left = self.radius() > start_edge.radius()
                        else:
                            left = True
                    if left:  # self is left
                        self.angle = np.pi * 2.0
                    else:
                        self.angle = 0.0
                elif self.is_line():
                    if start_lefthand:
                        self.angle = 0.0
                    else:
                        self.angle = np.pi * 2.0
            else:  # start is a line
                if self.is_arc():
                    if self.n1_direction_righthand():
                        self.angle = np.pi * 2.0
                    else:
                        self.angle = 0.0

        logger.debug("end set_direction_angle: angle is %s", self.angle)

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
        start_edge.log_edge("===> start")
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
                        arc_edge.set_direction_angle(start_edge)
                        left = not arc_edge.arc_line_direction_lefthand(start_edge, self, builder)
                        log_lefthand(left, self, arc_edge)
                        if left:  # self is left
                            myself_angle = 0.0
                        else:
                            myself_angle = np.pi * 2.0
                        left = myself_angle > other_angle

                        log_lefthand(left, self, nbr_edge)
                        logger.debug("#1a: end of myself_direction_lefthand: ==> %s", left)
                        return left
    
            elif self.is_arc():
                if not ignore_start:
                    if self.line_arc_close_together(start_edge):
                        logger.debug("START LINE and SELF ARC close")
                        line_edge = start_edge.get_reverse_edge()  # reverse start edge
                        line_edge.log_name("REVERSE")
                        line_edge.set_direction_angle(start_edge)
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

        if is_same_angle(0.0, other_angle, atol=0.008):  # 1/2 degree
            # 360 or 0 degrees => turn 180 degrees
            logger.debug("-- ATTENTION: other %s turns nearly 180 degrees", nbr_edge.classname())
            logger.debug("   the angle is %s", other_angle)
            if start_edge.is_arc():
                if nbr_edge.is_line() and not ignore_start:
                    if start_edge.line_arc_close_together(nbr_edge):
                        logger.debug("START ARC and SELF LINE close")
                        arc_edge = start_edge.get_reverse_edge()  # reverse start edge
                        arc_edge.set_direction_angle(start_edge)
                        left = arc_edge.arc_line_direction_lefthand(start_edge, nbr_edge, builder)
                        log_lefthand(left, self, nbr_edge)
                        logger.debug("#5: end of myself_direction_lefthand: ==> %s", not left)
                        return not left

            elif nbr_edge.is_arc():
                if not ignore_start:
                    if nbr_edge.line_arc_close_together(start_edge):
                        logger.debug("START LINE and NEIGHBOR ARC close")
                        line_edge = start_edge.get_reverse_edge()  # reverse start edge
                        line_edge.set_direction_angle(start_edge)
                        left = nbr_edge.arc_line_direction_lefthand(start_edge, line_edge, builder)
                        logger.debug("#6: end of myself_direction_lefthand: ==> %s", left)
                        return left

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
        if not self.is_arc():
            logger.critical("FATAL: unexpected %s at position %s", self.classname(), self.n1)
            sys.exit(1)

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
        next_edge.set_direction_angle(start_edge)
        left = self.myself_direction_lefthand(start_edge,
                                              next_edge,
                                              builder,
                                              ignore_start=True)
        logger.debug("end of arc_line_direction_lefthand: arc lefthand = %s", left)
        return left

    def get_alpha_points(self, prev_n):
        n1 = self.n1
        n2 = self.n2
        a = 0.0
        if self.is_arc():
            angle = self.element.get_angle_of_arc()
            if angle > np.pi * 0.75:
                n = self.element.center_of_connection()
                a = normalise_angle(alpha_points(prev_n,
                                                 n1,
                                                 n))
                n1 = n

        a += normalise_angle(alpha_points(prev_n,
                                          n1,
                                          n2))
        return a, n1

class AreaBuilder(object):
    def __init__(self,
                 geom=None,
                 nolog=True,
                 rtol=1e-04,
                 atol=1e-04,
                 ndec=6):
        assert(geom is not None)
        self.rtol = rtol
        self.atol = atol
        self.ndec = ndec
        self.geom = geom
        self.area_list = []
        self.errors = 0
        self.journal = getJournal()
        self.nolog = nolog
        self.num_edges = self.geom.number_of_edges()

    def __str__(self):
        return "rtol: {}\n".format(self.rtol) + \
               "atol: {}\n".format(self.atol)

    def all_neighbors(self, n):
        for n in self.geom.g.neighbors(n):
            yield n

    def set_edge_attributes(self):
        self.geom.set_edge_attributes()

    def append_new_area(self, area_list, elements):
        a = Area(elements, self.geom.center, 0.0)
        a.set_type(TYPE_AIR)  # air
        area_list.append(a)
        return a

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
        if self.nolog:
            disable_logging()
        self.set_edge_attributes()
        
        for n in self.geom.nodes():
            finished = False
            while not finished:
                finished = True
                nbrs = [nbr for nbr in self.geom.g.neighbors(n)]
                for next_n in nbrs:
                    result = self.get_new_area(n, next_n)
                    if result['ok']:
                        a = self.append_new_area(self.area_list, result['area'])
                        logger.debug("Area %s found", a.identifier())
        if self.nolog:
            enable_logging()
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
            result['circle'] = True
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

        prev_node = info_curr.n1
        next_n1 = info_next.n1
        next_n2 = info_next.n2

        alpha, prev_node = info_next.get_alpha_points(prev_node)

        logger.debug("--> alpha %s <--", alpha)
        c = 1
        while not (nodes_are_equal(next_n1, start_n1) and
                   nodes_are_equal(next_n2, start_n2)):
            c += 1
            if c > self.num_edges * 2:
                logger.error("FATAL: *** over %s elements in area ? ***",
                             self.num_edges)
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

            next_n1 = info_next.n1
            next_n2 = info_next.n2

            a, prev_node = info_next.get_alpha_points(prev_node)
            alpha += a
            logger.debug("--> alpha %s (+ %s) <--", alpha, a)

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

        if len(nbrs) == 1:
            logger.debug("end of next_edge_lefthand_side: one neighbor")
            return nbrs[0]

        logger.debug("-- NODE %s has %s neighbors", start_edge.n2, len(nbrs))
        logger.debug("-- startangle is %s", start_edge.startangle)

        all_nbrs = []
        logger.debug("candidates are:") 
        for i, nbr in enumerate(nbrs):
            nbr.set_direction_angle(start_edge)
            nbr.log_name("EDG-{}".format(i))
            nbr.log_edge("++++ candidate")
            all_nbrs.append((nbr.angle, i, nbr))
        logger.debug("--")
        all_nbrs.sort()

        a1, i1, nbr1 = all_nbrs[0]  # lefthandside
        nbr1.log_edge("==== first")
        for a2, i2, nbr2 in all_nbrs[1:]:
            nbr2.log_edge("==== next")
            if nbr2.myself_direction_lefthand(start_edge, nbr1, self):
                nbr1 = nbr2
                logger.debug("-- %s ist neuer lefthand neighbor", nbr1.name)
            else:
                logger.debug("-- %s bleibt lefthand neighbor", nbr1.name)
        nbr1.log_edge(">>>> LEFT")
        logger.debug("end of next_edge_lefthand_side")
        return nbr1

    def create_inner_corner_auxiliary_areas(self, startangle, endangle):
        logger.debug("begin of create_inner_corner_auxiliary_areas")
        if not self.geom.is_inner:
            logger.debug("end of create_inner_corner_auxiliary_areas: not inner")
            return False

        self.set_edge_attributes()

        start_nodes = [n for n in self.geom.angle_nodes(
            self.geom.center,
            startangle, 1e-3, 1e-3)]
        for n in start_nodes:
            nbrs = self.geom.get_neighbors(n)
            logger.debug("start node %s has %s neighbors", n, len(nbrs))
        logger.debug("corner nodes: %s", self.geom.start_corners)
        logger.debug("end nodes: %s", self.geom.end_corners)
        start_cp, start_exists = self.geom.get_start_airgap_corner()
        end_cp, end_exists = self.geom.get_end_airgap_corner()

        if start_exists and end_exists:
            logger.debug("end of create_inner_corner_auxiliary_areas: no aktion")
            return False

        airgap_line, airgap_el = self.get_inner_airgap_line()
        if not airgap_el:
            logger.debug("end of create_inner_corner_auxiliary_areas: no airgapline found")
            return False

        logger.debug("airgapline found !!")
        airgap_nodes = [n for n in airgap_line[1:]]
        del airgap_nodes[-1]

        created = False
        if not start_exists:
            cp = self.geom.start_corners[-1]
            logger.debug("Start Corner: %s -- %s", cp, start_cp)
            start_line = Line(Element(start=cp, end=start_cp),
                              color='red',
                              linestyle='dotted')

            start_cp = start_line.node2(self.ndec)
            for n in airgap_nodes:
                ag_line = Line(Element(start=start_cp, end=n))
                points = self.geom.get_intersection_points(airgap_el, ag_line, n)
                if not points:  # no intersection
                    d = distance(self.geom.center, n)
                    if np.isclose(d, self.geom.max_radius):
                        self.geom.add_arc(start_cp,
                                          n,
                                          self.geom.center,
                                          self.geom.max_radius,
                                          color='red',
                                          linestyle='dotted')
                    else:
                        self.geom.add_line(start_cp,
                                           n,
                                           color='red',
                                           linestyle='dotted')
                    created = True
                    start_node = self.geom.get_node(start_cp)
                    self.geom.add_edge(cp, start_node, start_line)
                    result = self.get_new_area(start_node, n)
                    if result['ok']:
                        self.append_new_area(self.geom.area_list,
                                             result['area'])
                    self.geom.set_start_corners(self.geom.center, 0.0)
                    break

        if not end_exists:
            cp = self.geom.end_corners[-1]
            logger.debug("End Corner: %s -- %s", cp, end_cp)
            end_line = Line(Element(start=cp, end=end_cp),
                            color='red',
                            linestyle='dotted')
            end_cp = end_line.node2(self.ndec)
            airgap_nodes.reverse()
            for n in airgap_nodes:
                ag_line = Line(Element(start=end_cp, end=n))
                points = self.geom.get_intersection_points(airgap_el, ag_line, n)
                if not points:  # no intersection
                    d = distance(self.geom.center, n)
                    if np.isclose(d, self.geom.max_radius):
                        self.geom.add_arc(n, end_cp,
                                          self.geom.center,
                                          self.geom.max_radius,
                                          color='red',
                                          linestyle='dotted')
                    else:
                        self.geom.add_line(end_cp, n,
                                           color='red',
                                           linestyle='dotted')
                    created = True
                    end_node = self.geom.get_node(end_cp)
                    self.geom.add_edge(cp, end_node, end_line)
                    result = self.get_new_area(n, end_node)
                    if result['ok']:
                        self.append_new_area(self.geom.area_list,
                                             result['area'])
                    self.geom.set_end_corners(self.geom.center, self.geom.alfa)
                    break

        logger.debug("end of create_inner_corner_auxiliary_areas")
        return created

    def get_airgap_line(self, start_node, end_node, area):
        logger.debug("get_airgap_line")

        self.set_edge_attributes()

        nodes = [n for n in area.list_of_nodes()]
        if not nodes:
            logger.debug("end of get_airgap_line: no nodes found")
            return [], []

        n1 = nodes[0]
        if points_are_close(start_node, n1):
            n2 = nodes[-1]
        else:
            n2 = n1
            for n1 in nodes[1:]:
                if points_are_close(start_node, n1):
                    break
                n2 = n1

        if not points_are_close(start_node, n1):
            logger.debug("end of get_airgap_line: not close to start-node")
            return [], []

        logger.debug("START EDGE FOUND: %s - %s", n1, n2)
        nodes = [n1, n2]
        info = self.get_edge_info(n1, n2)
        elements = [info.element]

        while not points_are_close(end_node, n2):
            info.set_start_angle()
            info = self.next_edge_lefthand_side(info)
            if not info:  # bad
                return [], []
            n2 = info.n2
            nodes.append(n2)
            elements.append(info.element)

        logger.debug("end of get_airgap_line #%s", len(nodes))
        return nodes, elements

    def get_inner_airgap_line(self):
        logger.debug("begin of get_inner_airgap_line")
        assert(self.geom.is_inner)
        assert(self.geom.area_list)

        area = [a for a in self.geom.area_list if a.close_to_ag_endcorner]
        if len(area) != 1:
            logger.debug("end of get_inner_airgap_line: %s areas found", len(area))
            return [], []

        start_node = self.geom.end_corners[-1]
        logger.debug("START NODE %s", start_node)
        end_node = self.geom.start_corners[-1]
        logger.debug("END NODE %s", end_node)

        return self.get_airgap_line(start_node, end_node, area[0])

    def close_outer_winding_areas(self):
        logger.debug("close_outer_winding_areas")

        airgap_line, airgap_el = self.get_outer_airgap_line()
        logger.debug("Outer Airgap with %s Nodes", len(airgap_line))

        if len(airgap_line) < 5:
            return False

        n1 = None
        dist_n1 = 0.0

        e_prev = None
        n_prev = airgap_line[0]
        dist_prev = distance(self.geom.center, n_prev)
        alpha_prev = alpha_line(self.geom.center, n_prev)
        alpha_start = alpha_prev

        lines_created = 0
        for n in airgap_line[1:]:
            dist = distance(self.geom.center, n)
            alpha = alpha_line(self.geom.center, n)
            if not n1:
                if dist > dist_prev and alpha < alpha_prev:
                    n1 = n_prev
                    dist_n1 = dist_prev
            else:
                if np.isclose(dist_n1, dist, rtol=1e-3, atol=1e-3):
                    line = Line(Element(start=n1, end=n))
                    if e_prev.intersect_line(line):
                        logger.debug("___LINE NOT POSSIBLE___")
                    else:
                        self.geom.add_line(n1, n, color='red')
                        lines_created += 1
                    n1 = None
                    dist_n1 = 0.0
            if not n1:
                e_prev = self.geom.get_edge_element(n_prev, n)
            n_prev = n
            dist_prev = dist
            alpha_prev = alpha

        return lines_created > 0

    def get_outer_airgap_line(self):
        logger.debug("begin of get_outer_airgap_line")
        assert(self.geom.is_outer)
        assert(self.geom.area_list)

        area = [a for a in self.geom.area_list if a.close_to_ag_startcorner]
        if len(area) != 1:
            logger.debug("end of get_outer_airgap_line: %s areas found", len(area))
            return [], []

        start_node = self.geom.start_corners[0]
        logger.debug("START NODE %s", start_node)
        end_node = self.geom.end_corners[0]
        logger.debug("END NODE %s", end_node)

        return self.get_airgap_line(start_node, end_node, area[0])

    def create_one_area_group(self, areas):
        logger.debug("begin of create_one_area_group")
        assert(len(self.area_list) == 0)

        self.set_edge_attributes()
        if not self.build_group_area(areas):
            logger.warning("Creation of an areagroup failed")
            return True
        return False

    def create_area_groups(self, list_of_areas):
        logger.debug("begin of create_area_groups")
        assert(len(self.area_list) == 0)

        area_id_list = [a.id for a in list_of_areas]
        timer = Timer(start_it=True)
        #   ---
        def delete_id(id):
            try:
                i = area_id_list.index(id)
            except ValueError:
                return True
            area_id_list[i] = 0
        #   ---

        group_list = {}
        for area in list_of_areas:
            if not area.id in area_id_list:
                continue

            areas = [area]
            delete_id(area.id)
            for a in list_of_areas:
                if area.id == a.id:
                    continue
                if not a.id in area_id_list:
                    continue
                if area.is_in_touch_with_area(self.geom, a):
                    areas.append(a)
                    delete_id(a.id)

            group_list[area.id] = areas

        group_keys = [int(x) for x in group_list.keys()]
        logger.debug("=== Area Group List ===")
        self.set_edge_attributes()
        errors = 0
        disable_logging()
        for k in group_keys:
            if not self.build_group_area(group_list[k]):
                logger.warning("Creation of an areagroup failed")
                errors += 1
        enable_logging()

        t = timer.stop("areagroups created in %0.4f seconds")
        logger.debug("end of create_area_groups: %s groups created",
                     len(self.area_list))
        return errors > 0

    def build_group_area(self, area_list):
        id_list = [a.id for a in area_list]
        logger.debug("Area Group: %s", id_list)

        max_x = 0
        area = None
        for a in area_list:
            if a.max_x > max_x:
                max_x = a.max_x
                area = a

        x0, y0 = -9999.0, 0.0
        for x, y in area.list_of_nodes():
            if x > x0:
                x0 = x
                y0 = y

        node0 = (x0, y0)
        if not self.geom.g.has_node(node0):
            node0 = self.geom.find_the_node(node0)
            x0, y0 = node0

        nbrs = self.geom.get_neighbors(node0)
        logger.debug("Neighbors of %s", node0)

        alpha0 = alpha_line(node0, (x, -9999.0))
        angle0 = np.pi * 2.0
        for n in nbrs:
            logger.debug("** nbr: %s", n)
            alpha1 = alpha_line(node0, n)
            angle1 = alpha_angle(alpha1, alpha0)
            logger.debug("** nbr=%s,  alpha=%s,  angle=%s", n, alpha1, angle1)
            logger.debug("** if %s < %s then", angle1, angle0)
            if angle1 < angle0:
                node1 = n
                angle0 = angle1

        logger.debug("** Start with nodes %s ==> %s", node0, node1)
        rslt = self.get_new_area(node0, node1)
        if rslt['ok'] and rslt.get('circle', False):
            rslt['ok'] = False
            rslt['reverse'] = True

        logger.debug("** Result: ok=%s, reverse=%s", rslt['ok'], rslt['reverse'])
        if not rslt['ok'] and rslt['reverse']:
            # new areagroup
            a = self.append_new_area(self.area_list, rslt['area'])
            a.areas_of_group = area_list
            return True
        return False
