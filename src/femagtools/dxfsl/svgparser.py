"""

Geometry Parser for SVG files

"""
import logging
import re
import lxml.etree as ET
from .shape import Circle, Arc, Line, Element
import numpy as np

logger = logging.getLogger(__name__)

def get_center(r, p1, p2, sweep):
    """return center point coordinates of arc"""
    dp = p2-p1
    s = np.linalg.norm(dp)
    delta = np.arctan2(dp[1], dp[0])
    if s < 2*r:
        if sweep == 0:
            alfa = delta - np.arctan2(np.sqrt(r**2-s**2/4), s/2)
        else:
            alfa = delta + np.arctan2(np.sqrt(r**2-s**2/4), s/2)
    else:
        alfa = delta
    return p1[0] + r*np.cos(alfa), p1[1] + r*np.sin(alfa)


def get_angles(sweep, center, p1, p2):
    x1, y1 = (p1-center)
    x2, y2 = (p2-center)
    if sweep == 0:
        return np.arctan2(y2, x2), np.arctan2(y1, x1)
    return np.arctan2(y1, x1), np.arctan2(y2, x2)


def get_shapes(path):
    """return list of node elements (A, L)"""
    state = ''
    p = []
    for s in [s for s in re.split('([AML])|,|\\s+',path) if s]:
        if state == '':
            state = s[0]
        elif state == 'M':
            p.append(float(s))
            if len(p) == 2:
                p1 = np.array(p)
                p = []
                state = ''
        elif state == 'L':
            p.append(float(s))
            if len(p) == 2:
                p2 = np.array(p)
                logger.debug("Line %s -> %s",
                            p1, p2)
                yield Line(Element(start=p1, end=p2))
                p1 = p2.copy()
                p = []
                state = ''
        elif state == 'A':
            p.append(float(s))
            if len(p) == 7:
                sweep = int(p[-3])
                p2 = np.array(p[-2:])
                r = p[0]
                center = get_center(r, p1, p2, sweep)
                start, end = get_angles(sweep, center, p1, p2)
                logger.debug("Arc center %s r %f %f -> %f",
                            center, r, start, end)
                yield Arc(Element(center=center,
                                  radius=r,
                                  start_angle=start*180/np.pi,
                                  end_angle=end*180/np.pi))
                p1 = p2.copy()
                p = []
                state = ''
        else:
            raise ValueError(f"unsupported path {state}")


def svgshapes(svgfile):
    svg = ET.parse(svgfile)
    for p in svg.findall(".//{http://www.w3.org/2000/svg}path"):
        yield from get_shapes(p.get('d'))
    for p in svg.findall(".//{http://www.w3.org/2000/svg}line"):
        yield Line(Element(start=[float(p.get('x1')), float(p.get('y1'))],
                           end=[float(p.get('x2')), float(p.get('y2'))]))
    for p in svg.findall(".//{http://www.w3.org/2000/svg}circle"):
        center = float(p.get('cx')), float(p.get('cy'))
        yield Circle(Element(center=center, radius=float(p.get('r'))))
