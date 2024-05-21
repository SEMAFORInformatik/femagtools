"""
 Geom Parser for femm files
"""

import numpy as np
from .shape import Arc, Line, Element
from .functions import distance, alpha_line, points_are_close, point_on_arc

fem_points = []


def read_fem_points(f, num):
    for x in range(num):
        p = f.readline().split()
        fem_points.append([float(p[0]), float(p[1])])


def read_fem_lines(f, num):
    for x in range(num):
        p = f.readline().split()
        i1 = int(p[0])
        i2 = int(p[1])
        p1 = fem_points[i1]
        p2 = fem_points[i2]
        if points_are_close(p1, p2):
            logger.warning("FEMM: Line with points close together")
            logger.warning("      p1 = %s, p2 =%s", p1, p2)
        yield Line(Element(start=p1, end=p2))


def read_fem_arcs(f, num):
    for x in range(num):
        p = f.readline().split()
        i1 = int(p[0])
        i2 = int(p[1])
        alpha = float(p[2])
        p1 = fem_points[i1]
        p2 = fem_points[i2]
        if points_are_close(p1, p2):
            logger.warning("FEMM: Arc with points close together")
            logger.warning("      p1 = %s, p2 = %s", p1, p2)
        for e in get_fem_arc(p1, p2, alpha):
            yield e


def get_fem_arc(pA, pB, alfa):
    alpha = alfa/180.0*np.pi/2.0
    y = distance(pA, pB) / 2.0
    x = y / np.tan(alpha)
    r = np.sqrt(x**2 + y**2)

    delta = alpha_line(pA, pB)

    c = [pA[0] + y, pA[1] + x]
    phi = alpha_line(pA, c) + delta
    pC = point_on_arc(pA, r, phi)

    startangle = alpha_line(pC, pA)
    endangle = alpha_line(pC, pB)
    yield Arc(Element(center=pC,
                      radius=r,
                      start_angle=startangle*180/np.pi,
                      end_angle=endangle*180/np.pi))


def femshapes(femfile):
    f = open(femfile, 'r')

    for data in f:
        text = data.split()
        if text[0] == '[NumPoints]':
            read_fem_points(f, int(text[2]))
        elif text[0] == '[NumSegments]':
            for e in read_fem_lines(f,  int(text[2])):
                yield e
        elif text[0] == '[NumArcSegments]':
            for e in read_fem_arcs(f,  int(text[2])):
                yield e
