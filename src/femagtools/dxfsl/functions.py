# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.functions
  ~~~~~~~~~~~~~~~~~~~~~~~~~~

  internal utility functions

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function
import logging
import numpy as np
import copy
import time
import multiprocessing

logger = logging.getLogger('femagtools.functions')


def less_equal(v1, v2, rtol=1e-3, atol=1e-8):
    if np.isclose(v1, v2, rtol, atol):
        return True
    return v1 < v2


def less(v1, v2, rtol=1e-3, atol=1e-8):
    if np.isclose(v1, v2, rtol, atol):
        return False
    return v1 < v2


def greater_equal(v1, v2, rtol=1e-3, atol=1e-8):
    if np.isclose(v1, v2, rtol, atol):
        return True
    return v1 > v2


def greater(v1, v2, rtol=1e-3, atol=1e-8):
    if np.isclose(v1, v2, rtol, atol):
        return False
    return v1 > v2


def alpha_line(center, p):
    return np.arctan2(p[1]-center[1], p[0]-center[0])


def alpha_points(a, b, c):
    alpha_a_b = alpha_line(a, b)
    alpha_b_c = alpha_line(b, c)
    return alpha_angle(alpha_a_b, alpha_b_c)


def alpha_triangle(a, b, c):
    if np.isclose(a, 0.0) or np.isclose(b, 0.0) or np.isclose(c, 0.0):
        return np.nan
    if a < 0.0 or b < 0.0 or c < 0.0:
        return np.nan
    if a + b < c:
        return np.nan
    if b + c < a:
        return np.nan
    if c + a < b:
        return np.nan

    cos_alpha = (a**2 - b**2 - c**2)/(-2*b*c)
    if cos_alpha >= 1.0:
        return np.nan
    if cos_alpha <= -1.0:
        return np.nan

    rslt = np.arccos(cos_alpha)
    if np.isnan(rslt):
        logger.debug("FATAL: arccos({}) yields nan.".format(cos_alpha))
        print('FATAL: arccos({}) yields nan. === {}/{}/{}'.
              format(cos_alpha, a, b, c))
    return rslt


def alpha_angle(startangle, endangle, rtol=1e-3, atol=1e-8):
    if less_equal(endangle, startangle, rtol=rtol, atol=atol):
        endangle += 2.0*np.pi
    angle = endangle - startangle
    if less_equal(angle, 2.0*np.pi):
        return angle
    return angle - 2.0*np.pi


def less_angle(alpha1, alpha2):
    angle = alpha_angle(alpha1, alpha2)
    return (angle < np.pi or angle > 2.0*np.pi)


def greater_angle(alpha1, alpha2):
    angle = alpha_angle(alpha1, alpha2)
    return not (angle < np.pi or angle > 2.0*np.pi)


def max_angle(alpha1, alpha2):
    if greater_angle(alpha1, alpha2):
        return alpha1
    return alpha2


def min_angle(alpha1, alpha2):
    if less_angle(alpha1, alpha2):
        return alpha1
    return alpha2


def point(center, radius, alpha, rnd=-1):
    if rnd >= 0:
        return (round(center[0]+radius*np.cos(alpha), rnd),
                round(center[1]+radius*np.sin(alpha), rnd))

    return (center[0]+radius*np.cos(alpha),
            center[1]+radius*np.sin(alpha))


def round_point(p, n):
    """round"""
    return (round(p[0], n), round(p[1], n))


def line_m(p1, p2, none_val=None, dec=0):
    if np.isclose(p2[0]-p1[0], 0.0):
        return none_val
    if dec > 0:
        return round((p2[1]-p1[1]) / (p2[0]-p1[0]), dec)
    return (p2[1]-p1[1]) / (p2[0]-p1[0])


def line_n(p, m):
    if m is None:
        return None
    return p[1] - m * p[0]


def distance(p1, p2):
    assert(len(p1) > 1)
    assert(len(p2) > 1)
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def area_size(points):
    if not points:
        return 0.0

    x = [p[0] for p in points]
    y = [p[1] for p in points]
    x.append(x[0])
    y.append(y[0])

    sz = 0.0
    for i in range(len(x) - 1):
        sz += x[i] * y[i+1] - y[i] * x[i+1]
    return np.absolute(sz) / 2


def middle_point_of_arc(center, radius, p1, p2, rtol=1e-3, atol=1e-8):
    alpha_p1 = alpha_line(center, p1)
    alpha_p2 = alpha_line(center, p2)
    if alpha_p1 < 0.0:
        alpha_p1 += np.pi * 2.0
    if alpha_p2 < 0.0:
        alpha_p2 += np.pi * 2.0

    if np.isclose(alpha_p1, alpha_p2, rtol, atol):
        return p1

    if greater_equal(alpha_p1, 0.0):
        if alpha_p2 < alpha_p1:
            alpha_p2 += 2.0*np.pi
    else:
        if less_equal(alpha_p2, alpha_p1):
            alpha_p1 += 2.0*np.pi

    if np.isclose(alpha_p1, alpha_p2):
        return copy.copy(p1)

    alpha_pm = (alpha_p1+alpha_p2) / 2.0
    return point(center, radius, alpha_pm)


def middle_angle(alpha1, alpha2):
    if alpha1 < 0.0:
        alpha1 += np.pi * 2.0
    if alpha2 < 0.0:
        alpha2 += np.pi * 2.0
    if alpha2 < alpha1:
        alpha2 += np.pi * 2.0
    return normalise_angle((alpha1 + alpha2) / 2.0)


def third_angle(alpha1, alpha2):
    a1 = normalise_angle(alpha1)
    a2 = normalise_angle(alpha2)

    if np.isclose(a1, a2):
        return a1

    if greater_equal(a1, 0.0):
        if a2 < a1:
            a2 += 2.0*np.pi
    else:
        if less_equal(a2, a1):
            a1 += 2.0*np.pi

    if np.isclose(a1, a2):
        return copy.copy(a1)

    return normalise_angle((a1+a2)/3.0)


def middle_point_of_line(p1, p2):
    return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2]


def lines_intersect_point(p_L1, m_L1, n_L1, p_L2, m_L2, n_L2):
    if m_L1 is None:
        return (p_L1[0], p_L2[1])
    if m_L2 is None:
        return (p_L2[0], p_L1[1])

    x = (n_L2-n_L1) / (m_L1-m_L2)
    y = m_L1 * x + n_L1
    return (x, y)


def intersect_point(p, L_p, L_m, L_n):
    if L_m is None:
        m = 0.0
    elif L_m == 0.0:
        m = None
    else:
        m = -1/L_m
    n = line_n(p, m)
    return lines_intersect_point(L_p, L_m, L_n, p, m, n)


def mirror_point(p, L_p, L_m, L_n):
    ps = intersect_point(p, L_p, L_m, L_n)
    x = p[0] - 2.0 * (p[0] - ps[0])
    y = p[1] - 2.0 * (p[1] - ps[1])
    return (x, y)


def points_are_close(p1, p2, rtol=1e-05, atol=1e-08):
    if not (p1 and p2):
        return False

    return (np.isclose(p1[0], p2[0], rtol, atol) and
            np.isclose(p1[1], p2[1], rtol, atol))


def point_greater_equal(p1, p2, rtol=1e-05, atol=1e-08):
    if not (p1 and p2):
        return False
    if np.isclose(p1[0], p2[0], rtol, atol):  # x equal
        return greater_equal(p1[1], p2[1], rtol, atol)  # y >= equal
    return greater_equal(p1[0], p2[0], rtol, atol)  # x >= equal


def nodes_are_equal(n1, n2):
    if not (n1 and n2):
        return False
    return (np.isclose(n1[0], n2[0], rtol=0.0, atol=1e-08) and
            np.isclose(n1[1], n2[1], rtol=0.0, atol=1e-08))


def point_in_region(p, x_min, x_max, y_min, y_max):
    if p[0] < x_min or p[0] > x_max:
        return False
    if p[1] < y_min or p[1] > y_max:
        return False
    return True


def within_interval(x, v1, v2, rtol=1e-4, atol=1e-8):
    """ returns true if x is in interval [v1, v2]
    """
    return np.logical_and(greater_equal(x, v1, rtol, atol),
                          less_equal(x, v2, rtol, atol))


def normalise_angle(alpha):
    """ returns angle alpha as value within interval [-pi, pi]
    """
    while alpha < -np.pi:
        alpha += 2*np.pi

    while alpha > np.pi:
        alpha -= 2*np.pi

    return alpha


def positive_angle(alpha):
    """ returns a positive value for angle alpha
    """
    while alpha < 0.0:
        alpha += 2*np.pi
    while alpha > 2*np.pi:
        alpha -= 2*np.pi
    return alpha


def elevation_angle(alpha):
    """ returns a positive angle for elevation
    """
    while alpha < 0.0:
        alpha += np.pi
    return alpha


def is_same_angle(angle1, angle2, atol=0.001):
    """ returns true if angles are equal
    """
    return (np.isclose(np.cos(angle1), np.cos(angle2), atol=atol) and
            np.isclose(np.sin(angle1), np.sin(angle2), atol=atol))


def part_of_circle(startangle, endangle, dec_place=2):
    """ returns the number of segments included in the circle
      if an integer number, 0 otherwise
    """
    start = normalise_angle(startangle)
    end = normalise_angle(endangle)

    if np.isclose(start, end):
        return 1

    if end > start:
        angle = end - start
    else:
        angle = 2*np.pi + end - start

    if angle != 0.0:
        x = float(round(2*np.pi/angle, dec_place))
    else:
        x = float(0.0)
    logger.debug("part_of_circle: {}".format(x))
    if x.is_integer():
        return int(x)
    return 0


def gcd(x, y):
    """ The function returns the greatest common divisor
    """
    while x > 0 and y > 0:
        if x >= y:
            x = x - y
        else:
            y = y - x
    return x+y


def is_angle_outside(startangle, endangle, alpha, rtol=1e-08, atol=1e-08):
    return not is_angle_inside(startangle, endangle, alpha, rtol=rtol, atol=atol)


def is_angle_inside(startangle, endangle, alpha, rtol=1e-08, atol=1e-08):
    start = normalise_angle(startangle)
    end = normalise_angle(endangle)
    mid = normalise_angle(alpha)

    if np.isclose(start, end, rtol=rtol, atol=atol):
        # In diesem Fall ist alles 'inside'
        return True
    if np.isclose(mid, start, rtol=rtol, atol=atol):
        return True
    if np.isclose(mid, end, rtol=rtol, atol=atol):
        return True

    if end < start:
        if mid > start:
            return True
        return mid < end
    else:
        if mid < start:
            return False
        return mid < end


def is_point_outside_region(p, center, inner_radius, outer_radius,
                            startangle, endangle,
                            rtol=1e-08, atol=1e-08):
    alpha = alpha_line(center, p)
    if is_angle_outside(startangle, endangle, alpha):
        return True
    dist = distance(center, p)
    return not within_interval(dist, inner_radius, outer_radius)


def is_point_inside_region(p, center,
                           inner_radius, outer_radius,
                           startangle, endangle,
                           rtol=1e-08, atol=1e-08):
    return not is_point_outside_region(p, center,
                                       inner_radius, outer_radius,
                                       startangle, endangle,
                                       rtol=rtol, atol=atol)


def get_angle_of_arc(startangle, endangle):
    if np.isclose(startangle, endangle, rtol=1e-06, atol=1e-09):
        return 2.0*np.pi  # is a circle

    start = normalise_angle(startangle)
    end = normalise_angle(endangle)

    if less(start, 0.0):
        start += 2.0*np.pi
    while less(end, start):
        end += 2.0*np.pi

    if np.isclose(start, end, rtol=1e-06, atol=1e-09):
        return 2.0*np.pi  # is a circle
    return end - start


def angles_on_arc(startangle, endangle, parts=8):
    alpha = get_angle_of_arc(startangle, endangle)
    num = max(int(alpha/(np.pi/parts)), 1)

    for x in range(0, num):
        yield float(x)/num*alpha + startangle
    if less(alpha, 2.0*np.pi):
        yield alpha + startangle


def point_on_arc(center, radius, angle):
    alpha = normalise_angle(angle)
    return (center[0] + radius * np.cos(alpha),
            center[1] + radius * np.sin(alpha))


def points_on_arc(center, radius, startangle, endangle, parts=8):
    start = normalise_angle(startangle)
    end = normalise_angle(endangle)
    for alpha in angles_on_arc(start, end, parts=parts):
        yield (center[0] + radius * np.cos(alpha),
               center[1] + radius * np.sin(alpha))


def points_on_line(p1, p2, parts=2):
    x_dist = (p2[0] - p1[0]) / parts
    y_dist = (p2[1] - p1[1]) / parts
    x, y = p1
    for i in range(1, parts):
        yield (x + i*x_dist, y + i*y_dist)


class Timer(object):
    def __init__(self, start_it=False):
        self.starttime = None
        if start_it:
            self.start()

    def __str__(self):
        if self.starttime is None:
            return "Timer is not running"
        stop = time.perf_counter()
        return "Timer is already running for %0.4f seconds" % (stop - self.starttime)

    def start(self):
        if self.starttime is not None:
            logger.error("Timer is already on")
        self.starttime = time.perf_counter()

    def stop(self, fmt=None, info=False):
        if self.starttime is None:
            logger.error("Timer is not running")
        stop = time.perf_counter()
        sec = stop - self.starttime
        self.starttime = None
        if fmt:
            if info:
                logger.info(fmt, sec)
            else:
                logger.debug(fmt, sec)
        return sec


class SimpleProcess(multiprocessing.Process):
    def __init__(self, name=None, no_processing=False):
        super(SimpleProcess, self).__init__()
        self._no_processing = no_processing
        self.name = name

    def without_processing(self):
        return self._no_processing

    def wait(self):
        if not self._no_processing:
            logger.debug("%s Waiting for Process", self.name)
            self.join()

    def start_task(self):
        if self._no_processing:
            self.run()
            return
        self.start()
