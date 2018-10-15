# -*- coding: utf-8 -*-
"""
    NOTE: This code is in highly experimental state.
          Use at your own risk.

  Author: Ronald Tanner
    Date: 2017/07/06
"""
from __future__ import print_function
import logging
import numpy as np
import copy

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


def alpha_angle(startangle, endangle):
    if less_equal(endangle, startangle):
        endangle += 2.0*np.pi
    angle = endangle - startangle
    if less_equal(angle, 2.0*np.pi):
        return angle
    return angle - 2.0*np.pi


def max_angle(alpha1, alpha2):
    angle = alpha_angle(alpha1, alpha2)
    if angle < np.pi or angle > 2.0*np.pi:
        return alpha2
    return alpha1


def min_angle(alpha1, alpha2):
    angle = alpha_angle(alpha1, alpha2)
    if angle < np.pi or angle > 2.0*np.pi:
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


def line_m(p1, p2):
    if np.isclose(p2[0]-p1[0], 0.0):
        return None
    return (p2[1]-p1[1]) / (p2[0]-p1[0])


def line_n(p, m):
    if m is None:
        return None
    return p[1] - m * p[0]


def distance(p1, p2):
    assert(len(p1) > 1)
    assert(len(p2) > 1)
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def middle_point_of_arc(center, radius, p1, p2, rtol=1e-3, atol=1e-8):
    alpha_p1 = alpha_line(center, p1)
    alpha_p2 = alpha_line(center, p2)

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

    return (a1+a2)/2.0


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

    return (a1+a2)/3.0


def middle_point_of_line(p1, p2):
    return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2]


def lines_intersect_point(p_L1, m_L1, n_L1, p_L2, m_L2, n_L2):
    if m_L1 is None:
        return (p_L1[0], p_L2[1])
    if m_L2 is None:
        return (p_L2[0], p_L1[1])

    x = (n_L2-n_L1) / (m_L1-m_L2)
    y = m_L1 * x + n_L1
    return [x, y]


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


def point_in_region(p, x_min, x_max, y_min, y_max):
    if p[0] < x_min or p[0] > x_max:
        return False
    if p[1] < y_min or p[1] > y_max:
        return False
    return True


def within_interval(x, v1, v2, rtol=1e-3, atol=1e-8):
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


def is_same_angle(angle1, angle2):
    """ returns true if angles are equal
    """
    return (np.isclose(np.cos(angle1), np.cos(angle2)) and
            np.isclose(np.sin(angle1), np.sin(angle2)))


def part_of_circle(startangle, endangle, pos=3):
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
        x = float(round(2*np.pi/angle, pos))
    else:
        x = float(0.0)
    logger.debug("part_of_circle: {}".format(x))
    if x.is_integer():
        return x
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


def is_angle_outside(startangle, endangle, alpha):
    return not is_angle_inside(startangle, endangle, alpha)


def is_angle_inside(startangle, endangle, alpha):
    start = normalise_angle(startangle)
    end = normalise_angle(endangle)
    mid = normalise_angle(alpha)

    if np.isclose(start, end, 1e-08):
        # In diesem Fall ist alles 'inside'
        return True
    if np.isclose(mid, start, 1e-08):
        return True
    if np.isclose(mid, end, 1e-08):
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
                            startangle, endangle):
    alpha = alpha_line(center, p)
    if is_angle_outside(startangle, endangle, alpha):
        return True
    dist = distance(center, p)
    return not within_interval(dist, inner_radius, outer_radius)


def is_point_inside_region(p, center,
                           inner_radius, outer_radius,
                           startangle, endangle):
    return not is_point_outside_region(p, center,
                                       inner_radius, outer_radius,
                                       startangle, endangle)


def get_angle_of_arc(startangle, endangle):
    if np.isclose(startangle, endangle):
        return 2.0*np.pi  # is a circle

    start = normalise_angle(startangle)
    end = normalise_angle(endangle)

    if less(start, 0.0):
        start += 2.0*np.pi
    while less(end, start):
        end += 2.0*np.pi

    if np.isclose(start, end):
        return 2.0*np.pi  # is a circle
    return end - start


def angles_on_arc(startangle, endangle, parts=8):
    alpha = get_angle_of_arc(startangle, endangle)
    num = max(int(alpha/(np.pi/parts)), 1)

    for x in range(0, num):
        yield float(x)/num*alpha + startangle
    if less(alpha, 2.0*np.pi):
        yield alpha + startangle


def points_on_arc(center, radius, startangle, endangle, parts=8):
    start = normalise_angle(startangle)
    end = normalise_angle(endangle)
    for alpha in angles_on_arc(start, end, parts=parts):
        yield (center[0] + radius * np.cos(alpha),
               center[1] + radius * np.sin(alpha))
