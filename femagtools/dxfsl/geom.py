"""manage a geometry with
    lines, circles and arcs built from DXF
    
  NOTE: This code is in highly experimental state.
        Use at your own risk.

  Author: Ronald Tanner
    Date: 2017/07/06
"""
import ezdxf
import dxfgrabber
import numpy as np
import numpy.linalg as la
import networkx as nx
import copy
import logging

logger = logging.getLogger(__name__)

ndec = 4  # number of decimals to round to

def less_equal(v1, v2, pickdist=1e-3):
    if np.isclose(v1, v2, pickdist):
        return True
    return v1 < v2

def greater_equal(v1, v2, pickdist=1e-3):
    if np.isclose(v1, v2, pickdist):
        return True
    return v1 > v2
    
def alpha_line(center, p):
    return np.arctan2(p[1]-center[1], p[0]-center[0])

def alpha_triangle(a, b, c):
    if np.isclose(a, 0.0) or np.isclose(b, 0.0) or np.isclose(c, 0.0):
        return float('nan')
    cos_alpha = (a**2 - b**2 - c**2)/(-2*b*c)
    if np.isnan(cos_alpha):
        return cos_alpha
    if a + b < c:
        return float('nan')
    return np.arccos(cos_alpha)

def point(center, radius, alpha):
    return (center[0]+radius*np.cos(alpha),
            center[1]+radius*np.sin(alpha))

def middle_of_arc(center, radius, p1, p2, pickdist=1e-3):
    alpha_p1 = alpha_line(center, p1)
    alpha_p2 = alpha_line(center, p2)
    if np.isclose(alpha_p1, alpha_p2, pickdist):
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

def middle_of_line(p1, p2):
    return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2]
    
def distance(p1, p2):
    assert(len(p1)>1)
    assert(len(p2)>1)
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def line_m(p1, p2):
    if np.isclose(p2[0]-p1[0], 0.0):
        return None
    return (p2[1]-p1[1]) / (p2[0]-p1[0])

def line_n(p, m):
    if m == None:
        return None
    return p[1] - m * p[0]

def lines_intersect_point(m_L1, n_L1, m_L2, n_L2):
    x = (n_L2-n_L1) / (m_L1-m_L2)
    y = m_L1 * x + n_L1
    return [x, y]
    
def points_are_close(p1, p2, rtol=1e-05, atol=1e-08):
    return np.isclose(p1[0], p2[0], rtol, atol) and \
           np.isclose(p1[1], p2[1], rtol, atol)

def in_range(x, v1, v2):
    """ Die Funktion prüft, ob der Wert x zwischen v1 und v1 liegt.
    """
    if not greater_equal(x, v1):
        return False
    if not less_equal(x, v2):
        return False
    return True

def normalise_angle(alpha):
    """ Die Funktion liefert den Winkel alpha als Wert zwischen - und + pi.
    """
    while alpha < -np.pi:
        alpha += 2*np.pi

    while alpha > np.pi:
        alpha -= 2*np.pi
        
    return alpha

def is_same_angle(angle1, angle2):
    """ Die Funktion prüft, ob die beiden Winkelmasse logisch gleich sind.
    """
    return np.isclose(np.cos(angle1), np.cos(angle2)) and \
           np.isclose(np.sin(angle1), np.sin(angle2))

def alpha_angle(startangle, endangle):
    if less_equal(endangle, startangle):
        endangle += 2.0*np.pi
    angle = endangle - startangle
    if less_equal(angle, 2.0*np.pi)         :
        return angle
    return angle - 2.0*np.pi

def part_of_circle(startangle, endangle):
    """ Die Funktion prüft, ob der Winkel ein ganzzahliger Teil eines
        Kreises ist und liefert den Nenner von 1/n.
    """
    start = normalise_angle(startangle)
    end = normalise_angle(endangle)
    
    if np.isclose(start, end):
        return 0
        
    if end > start:
        angle = end - start
    else:
        angle = 2*np.pi + end - start
    
    if angle != 0.0:
        x = float(round(2*np.pi/angle, 3))
    else:
        x = float(0.0)
    print("Teiler = {}".format(x))
    if x.is_integer():
        return x
    return 0
    
def is_outside(startangle, endangle, alpha):
    return not is_inside(startangle, endangle, alpha)

def is_inside(startangle, endangle, alpha):
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

def is_outside_region(p, center, inner_radius, outer_radius, startangle, endangle):
    alpha = alpha_line(center, p)
    if is_outside(startangle, endangle, alpha):
        return True
    dist = distance(center, p)
    return not in_range(dist, inner_radius, outer_radius)
    
def is_inside_region(p, center, inner_radius, outer_radius, startangle, endangle):
    return not is_outside_region(p, center, inner_radius, outer_radius, startangle, endangle)

def angles_on_arc(startangle, endangle):
    circle = np.isclose(startangle, endangle)
    if circle:
        endangle += 2.0*np.pi
    elif greater_equal(startangle, 0.0):
        if endangle < startangle:
            endangle += 2.0*np.pi
    else:
        if less_equal(endangle, startangle):
            startangle += 2.0*np.pi

    alpha = endangle - startangle

    num = int(alpha/(np.pi/8))
    
    for x in range(0, num):
        yield x/num*alpha + startangle
    if not circle:
        yield alpha + startangle
        
def points_on_arc(center, radius, startangle, endangle):
    start = normalise_angle(startangle)
    end = normalise_angle(endangle)
    for alpha in angles_on_arc(start, end):
        yield (center[0] + radius * np.cos(alpha), center[1] + radius * np.sin(alpha))

def reshape(geom, incl_bnd=False):
    """take a symmetry section"""
    n = round(np.pi/geom.argmax())
    logger.info("n %d", n)
    alpha = np.pi/n/2
    sector = geom.copy(0.0, alpha, incl_bnd=incl_bnd)
    if not sector.g:
        return sector
    
    # remove boundary edges
    edges = []
    d = sorted((sector.diameters[0], sector.diameters[-1]))
    for e in sector.g.edges():
        r = np.linalg.norm(e, axis=1)
        if ((np.isclose(e[0][1], 0, 1e-3) and np.isclose(e[1][1], 0, 1e-3)) or
            np.isclose(r, d[0]/2, 1e-3).all() or
            np.isclose(r, d[1]/2, 1e-3).all()):
            edges.append(e)
    for e in edges:
        sector.g.remove_edge(e[0], e[1])

    # add boundary arcs
    for dx in d:
        arc = Arc(Element(
            start_angle=0.0, end_angle=alpha*180/np.pi,
            center=(0.0, 0.0), radius=dx/2))
        a = [np.arctan2(p[1], p[0])
             for p in convex_hull(sector.g.nodes())
             if (np.isclose(np.linalg.norm(p), arc.radius) and
                 arc.is_inside(np.arctan2(p[1], p[0]), 1e-3))]
        a.insert(0, arc.startangle)
        a.append(arc.endangle)
        for i in range(len(a[:-1])):
            sector.add_arc(
                start_angle=a[i],
                end_angle=a[i+1],
                center=(0.0, 0.0), radius=dx/2)
            
    sector.add_line((d[0]/2, 0.0), (d[1]/2, 0.0))
    rot = np.array([[np.cos(alpha), -np.sin(alpha)],
                    [np.sin(alpha), np.cos(alpha)]])
    p0, p1 = np.dot(rot, (d[0]/2, 0)), np.dot(rot, (d[1]/2, 0))
    sector.add_line(p0, p1)

    return sector

def find_corners(nodes, all=False):
    """find corners of nodes"""
    if nodes:
        a = np.asarray(nodes).T
        if all:
            args = np.arctan2(a[1], a[0])
            rads2 = np.sum(a**2, axis=0)

            phimin, phimax = np.amin(args), np.amax(args)
            r2min, r2max = np.amin(rads2), np.amax(rads2)

            rindx = np.union1d(np.isclose(rads2, r2min, 1e-2).nonzero()[0],
                               np.isclose(rads2, r2max, 1e-2).nonzero()[0])
            phindx = np.union1d(np.isclose(args, phimin).nonzero()[0],
                                np.isclose(args, phimax).nonzero()[0])
            return np.asarray(nodes)[np.union1d(rindx, phindx)].tolist()
        else:
            # collect real corners only (TODO: there should be some simpler way)
            corners = []
            a = np.asarray(nodes)
            minX, minY = np.amin(a, axis=0)
            maxX, maxY = np.amax(a, axis=0)
            x = a[np.isclose(a.T[0], minX).nonzero()[0]]
            s = np.lexsort(x.T)
            corners.append(x[s[0]])
            if len(s) > 1:
                corners.append(x[s[-1]])
            x = a[np.isclose(a.T[0], maxX).nonzero()[0]]
            s = np.lexsort(x.T)
            corners.append(x[s[0]])
            if len(s) > 1:
                corners.append(x[s[-1]])

            y = a[np.isclose(a.T[1], minY).nonzero()[0]]
            s = np.lexsort((y.T[1], y.T[0]))
            corners.append(y[s[0]])
            if len(s) > 1:
                corners.append(y[s[-1]])
            y = a[np.isclose(a.T[1], maxY).nonzero()[0]]
            s = np.lexsort((y.T[1], y.T[0]))
            corners.append(y[s[0]])
            if len(s) > 1:
                corners.append(y[s[-1]])
            return set([tuple(c) for c in corners])

    return []


def remove_corners(self, g):
    """removes the corner nodes with their edges"""
    corners = [n for n in find_corners(g.nodes()) if g.degree(n) < 3]
    logger.debug("removing corners %s", corners)
    g.remove_nodes_from(corners)

def add_and_split(inp_elements, pickdist):
    out_elements = []

    for inp_el in inp_elements:
        new_elements = add_and_split_element(inp_el, out_elements, pickdist)
        out_elements += new_elements
    return out_elements

def get_elements(elements):
    for el in elements:
        if isinstance(el, list):
            yield get_elements(el)
        else:
            yield el
    
def add_and_split_element(el, elements, pickdist):
    for x in range(len(elements)):
        e = elements[x]
        assert(isinstance(e, Shape))
            
    return [el]
    
def add_or_join(g, n1, n2, entity, pickdist):
    """ adds a new entity to graph or joins entity with existing
    g: graph
    n1, n2: nodes
    entity
    """
    if isinstance(entity, Arc):
        for e1, e2, a in g.edges([n1, n2], data=True):
            if (e1 == n1 or e1 == n2) and \
               (e2 == n1 or e2 == n2):
                o = a['object']
                if isinstance(o, Arc) and \
                   np.isclose(o.center,
                              entity.center, pickdist).all() and \
                    np.isclose(o.radius,
                               entity.radius, pickdist) and \
                    not np.isclose(o.startangle, entity.startangle) and \
                    not np.isclose(o.endangle, entity.endangle):
                    logger.info("**** Circle %s %g (%g, %g) --> (%g, %g)",
                                o.center, o.radius, o.startangle, o.endangle,
                                entity.startangle, entity.endangle)
                    g.remove_edge(e1, e2)
                    g.add_node(o.center,
                               object=Circle(Element(center=o.center,
                                                     radius=o.radius)))
                    return

    g.add_edge(n1, n2, object=entity)
    

def get_nodes_of_paths(g, c):
    nodes = set()
    i = 1
    for s in c[:-1]:
        for e in c[i:]:
            logger.info("%s --%s:", s, e)
            for p in nx.all_simple_paths(g, s, e):
                logger.info("   %s", len(list(p)))
                nodes |= set([n for n in p if g.degree(n) < 3])
        i += 1
    return nodes


def convex_hull(nodes):
    """collect edges that build the convex hull."""

    # Get a local list copy of the nodes and sort them lexically.
    points = list(nodes)
    points.sort()
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of
    # their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        dx = a[0] - o[0], a[1] - o[1]
        dy = b[0] - o[0], b[1] - o[1]
        return dx[0] * dy[1] - dx[1] * dy[0]

    # Build lower hull
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the
    # beginning of the other list.
    return lower[:-1] + upper[:-1]


def remove_hull(g, pickdist):
    hull = sorted(convex_hull(g.nodes()))
    rmax = la.norm(hull[-1])
    outarc = [h for h in hull
              if np.isclose(la.norm(h),
                            rmax, pickdist)]
    logger.info("remove hull outarc %d", len(outarc))
    g.remove_nodes_from(outarc)
    for o in outarc:
        hull.remove(o)
    argmax = np.max([np.arctan2(h[1], h[0]) for h in hull])
    secvec = np.array((np.cos(argmax),
                       np.sin(argmax)))
    outsec = [h for h in g.edges() if np.isclose(
        h[0],
        la.norm(h[0])*secvec,
        pickdist).all() and np.isclose(
            h[1],
            la.norm(h[1])*secvec,
            pickdist).all()
    ]
    if outsec:
        logger.info("remove hull outsec %d (argmax %f)",
                    len(outsec), argmax)
        g.remove_edges_from(outsec)

    lostnodes = [n for n in g.nodes() if g.degree(n) == 0]
    logger.info("lostnodes %d", len(lostnodes))
    if lostnodes:
        g.remove_nodes_from(lostnodes)


def ccw(a, b, c):
    """
    checks if 3 points are placed counter clockwise or in line
    """
    cp = np.cross(np.asarray(b)-np.asarray(a),
                  np.asarray(c)-np.asarray(a))
    if np.allclose(cp, (0.0,)):
        return 0
    if cp > 0:
        return 1
    return -1


def intersect(a, b):
    """checks if 2 lines intersect"""
    return ccw(a.p1, a.p2, b.p1) != ccw(a.p1, a.p2, b.p2) and \
        ccw(b.p1, b.p2, a.p1) != ccw(b.p1, b.p2, a.p2)


class Element(object):
    """value object class"""
    def __init__(self, **kwargs):
        for k in kwargs.keys():
            setattr(self, k, kwargs[k])


def polylines(entity):
    """returns a collection of bulged vertices
    http://www.afralisp.net/archive/lisp/Bulges1.htm
    """
    if isinstance(entity.points, list):
        points = [(p[0], p[1]) for p in entity.points]
    else:
        points = [(p[0], p[1]) for p in entity.points()]
    i = 0
    if isinstance(entity.vertices, list):
        vertices = entity.vertices
    else:
        vertices = entity.vertices()
        
    for v in vertices:
        if hasattr(v, 'bulge'):
            b = v.bulge
        else:
            b = v.get_dxf_attrib('bulge', 0.0)
        p1 = points[i]
        try:
            p2 = points[i+1]
        except:
            p2 = points[0]
        if b != 0.0:
            dx, dy = p2[0] - p1[0], p2[1] - p1[1]
            c = np.sqrt(dx**2 + dy**2)
            s = b * c/2
            r = ((c/2)**2 + s**2)/2/s
            g = np.arctan2(dx, dy)
            a = 2*np.arctan(b)
            pc = (p1[0] + r*np.sin(g - np.pi/2 + a),
                  p1[1] + r*np.cos(g - np.pi/2 + a))
            if b < 0:
                sa = np.arctan2(p2[1]-pc[1], p2[0]-pc[0])
                ea = np.arctan2(p1[1]-pc[1], p1[0]-pc[0])
            else:
                sa = np.arctan2(p1[1]-pc[1], p1[0]-pc[0])
                ea = np.arctan2(p2[1]-pc[1], p2[0]-pc[0])
            logger.debug("Poly p1 %s p2 %s r %s", p1, p2, r)
            yield Arc(Element(center=pc, radius=np.abs(r),
                              start_angle=sa/np.pi*180,
                              end_angle=ea/np.pi*180))
        else:
            yield Line(Element(start=p1, end=p2))
        i += 1
     
    
def dxfshapes0(dxffile):
    """returns a collection of dxf entities (ezdxf)"""
    dwg = ezdxf.readfile(dxffile)
    id = 0
    # $ACADVER: AC1006 = R10, AC1009 = R11 and R12, AC1012 = R13, AC1014 = R14 AC1015 = Release 2000/0i/2
    # check units:
    # dwg.header['$ANGDIR'] 1 = Clockwise angles, 0 = Counterclockwise
    # dwg.header['$AUNITS'] Decimal Degrees, Deg/Min/Sec, Grads, Radians
    # dwg.header['$INSUNIT'] 1 = Inches; 2 = Feet; 3 = Miles; 4 = Millimeters; 5 = Centimeters; 6 = Meters
    # dwg.header['$LUNITS']
    for e in dwg.modelspace():
        if e.dxftype() == 'ARC':
            yield Arc(e.dxf)
        elif e.dxftype() == 'CIRCLE':
            logger.info("C %s, R %f", e.center[:2], e.radius)
            yield Circle(e.dxf)
        elif e.dxftype() == 'LINE':
            yield Line(e.dxf)
        elif e.dxftype() == 'POLYLINE':
            for p in polylines(e):
                yield p
        id += 1


def dxfshapes(dxffile, layers=[]):
    """returns a collection of dxf entities (dxfgrabber)"""
    dwg = dxfgrabber.readfile(dxffile)
    id = 0
    # $ACADVER: AC1006 = R10, AC1009 = R11 and R12, AC1012 = R13, AC1014 = R14 AC1015 = Release 2000/0i/2
    # check units:
    # dwg.header['$ANGDIR'] 1 = Clockwise angles, 0 = Counterclockwise
    # dwg.header['$AUNITS'] Decimal Degrees, Deg/Min/Sec, Grads, Radians
    # dwg.header['$INSUNIT'] 1 = Inches; 2 = Feet; 3 = Miles; 4 = Millimeters; 5 = Centimeters; 6 = Meters
    # dwg.header['$LUNITS']
    for e in dwg.modelspace():
        if not layers or e.layer in layers:
            if e.dxftype == 'ARC':
                yield Arc(e)
            elif e.dxftype == 'CIRCLE':
                logger.info("C %s, R %f", e.center[:2], e.radius)
                yield Circle(e)
            elif e.dxftype == 'LINE':
                yield Line(e)
            elif e.dxftype == 'POLYLINE':
                for p in polylines(e):
                    yield p
            id += 1

#############################
#           Corner          #
#############################

class Corner(object):
    def __init__(self, center, p):
        self.__p = p
        self.__dist = distance(center, p)
        self.__keep = False
       
    def point(self):
        return self.__p
        
    def is_equal(self, p):
        return np.isclose(self.__p[0], p[0]) and np.isclose(self.__p[1], p[1])
        
    def set_keep_node(self):
        self.__keep = True

    def keep_node(self):
        return self.__keep

    def __lt__(self, c):               
        if self.__dist < c.__dist:
            return 1
        else:
            return 0
        
    def __str__(self):
        return "Corner: p={}".format(self.__p)

#############################
#       Shape (Basis)       #
#############################

class Shape(object):
    """an abstract geometry with 2 points"""
    def start(self):
        return self.p1
    
    def end(self):
        return self.p2

    def xmin(self):
        return min(self.p1[0], self.p2[0])
    
    def xmax(self):
        return max(self.p1[0], self.p2[0])
    
    def ymin(self):
        return min(self.p1[1], self.p2[1])
    
    def ymax(self):
        return max(self.p1[1], self.p2[1])

    def dx(self):
        return (self.p2[0]-self.p1[0])
    
    def dy(self):
        return (self.p2[1]-self.p1[1])

    def m(self):
        return line_m(self.p1, self.p2)

    def n(self, m):
        return line_n(self.p1, m)
    
    def move(self, dist):
        self.p1 = self.p1[0] + dist[0], self.p1[1] + dist[1]
        self.p2 = self.p2[0] + dist[0], self.p2[1] + dist[1]

    def scale(self, factor):
        self.p1 = factor*self.p1[0], factor*self.p1[1]
        self.p2 = factor*self.p2[0], factor*self.p2[1]
        
    def transform(self, T, **kwargs):
        n = T.dot(np.array((self.p1[0], self.p1[1])))
        self.p1 = (n[0], n[1])
        n = T.dot(np.array((self.p2[0], self.p2[1])))
        self.p2 = (n[0], n[1])
        return self

    def split_shape(self, e, pickdist):
        return [self], [e]

    def __str__(self):
        return " {}/{}".format(self.p1, self.p2)
    
#############################
#       Circle (Shape)      #
#############################

class Circle(Shape):
    """a circle with center and radius"""
    def __init__(self, e):
        self.center = e.center[:2]
        self.radius = e.radius
        self.p1 = self.center[0]-self.radius, self.center[1]
        self.p2 = self.center[0]+self.radius, self.center[1]
        
    def render(self, renderer, color='blue', with_nodes=False):
        renderer.circle(self.center, self.radius, color)
        if with_nodes:
            renderer.point(self.center, 'ro', 'white')

    def move(self, dist):
        super(Circle, self).move(dist)
        self.center = self.center[0]+dist[0], self.center[1]+dist[1]

    def minmax(self):
        """ Die Funktion bestimmt das Minimum und Maximum auf der x- und der
            y-Achse (return [<min-x>, <max-x>, <min-y>, <max-y>])
        """
        return [self.center[0]-self.radius,self.center[0]+self.radius,
                self.center[1]-self.radius,self.center[1]+self.radius]

    def minmax_from_center(self, center):
        """ Die Funktion ermittelt den minimalen und maximalen Abstand vom Center 
        """
        d = distance(center, self.center)
        return (max(0.0, d - self.radius), d + self.radius)
        
    def get_nodes(self):
        """ Die Funktion liefert eine Liste von virtuellen Nodes, welche man
            zum Rechnen der convex_hull() benötigt.
        """
        return (p for p in points_on_arc(self.center, self.radius, 0.0, 0.0))
        
    def scale(self, factor):
        super(Circle, self).scale(factor)
        self.center = factor*self.center[0], factor*self.center[1]
        self.radius = factor*self.radius

    def transform(self, T, **kwargs):
        super(Circle, self).transform(T)
        n = T.dot(np.array((self.center[0], self.center[1])))
        self.center = (n[0], n[1])
        return self

    def center_of_connection(self):
        return (self.center[0] + self.radius, self.center[1])

    def intersect(self, line, pickdist=1e-3):
        """calculate intersection of this circle with line"""
        logger.debug(self)
        if np.isclose(line.dx(), 0):  # this is a vertical line
            x = line.p1[0]
            d2 = self.radius**2 - (x-self.center[0])**2
            if d2 >= 0:
                d = np.sqrt(d2)
                if np.isclose(d, 0):
                    return ((x, self.center[1]),)
                return ((x, d + self.center[1]),
                        (x, -d + self.center[1]))
            return []

        m = line.dy()/line.dx()
        q = line.p1[1]-self.center[1] - m*(line.p1[0] - self.center[0])
        A = 1 + m**2
        B = 2*m*q
        C = q**2 - self.radius**2

        d2 = B**2 - 4*A*C
        if d2 > 0:
            d = np.sqrt(d2)
            x = (-B + d)/2/A, (-B - d)/2/A
            if np.isclose(d, 0):
                return ((x[0]+self.center[0],
                         m*x[0] + q + self.center[1]),)
            return ((x[0] + self.center[0],
                     m*x[0] + q + self.center[1]),
                    (x[1] + self.center[0],
                     m*x[1] + q + self.center[1]))
        return []

    def intersect_line(self, line, pickdist=1e-3, include_end=False):
        """ Von einem Circle-Objekt und einem Line-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben.
        """
        m_L = line.m()
        m_c = 0.0
        p = []
        if m_L == None:
            p = [line.p1[0], self.center[1]]
        elif np.isclose(m_L, 0.0, pickdist):
            p = [self.center[0], line.p1[1]]
        else:
            m_c = -1/m_L
            p = lines_intersect_point(m_L, line.n(m_L), m_c, line_n(self.center, m_c))

        d = distance(self.center, p)
            
        if np.isclose(d, self.radius, pickdist):
            # Wenn der Abstand d dem Radius entspricht, handelt es sich um
            # eine Tangente und es gibt genau einen Schnittpunkt
            if include_end:
                return [p]
            else:
                return []
        if self.radius < d:
            # d liegt ausserhalb des Kreises -> kein Schnittpunkt
            return []
            
        A = np.sqrt(self.radius**2 - d**2)        
        delta = alpha_line(line.p1, line.p2)
        
        p1 = point(p, -A, delta)
        p2 = point(p, A, delta)

        # Die Schnittpunkte p1 und p2 sind bestimmt. Nun muss noch sicher
        # gestellt werden, dass sie innerhalb des Start- und Endpunkts der
        # Linie liegen
        if line.is_inside(p1, pickdist, include_end):
            if line.is_inside(p2, pickdist, include_end):
                return [p1, p2]
            else:
                return[p1]
        else:
            if line.is_inside(p2, pickdist, include_end):
                return[p2]
            else:
                return []
                
    def intersect_circle(self, circle, pickdist=1e-3, include_end=False):
        """ Von zwei Circle-Objekten werden die Schnittpunkte bestimmt
            und in einer Liste ausgegeben
        """
        d = distance(self.center, circle.center)
        arc = alpha_triangle(circle.radius, self.radius, d)
        if np.isnan(arc):
            return []
        arc_C = alpha_line(self.center, circle.center)
        p1 = point(self.center, self.radius, arc_C+arc)
        p2 = point(self.center, self.radius, arc_C-arc)
        if points_are_close(p1, p2, pickdist):
            # Tangente
            if include_end:
                return [p1]
            else:
                return []
        return [p1, p2]
        
    def intersect_arc(self, arc, pickdist=1e-3, include_end=False):
        """ Von einem Circle-Objekt und einem Arc-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben
        """
        assert(isinstance(arc, Arc))
        # Die Arbeit übernimmt das Arc-Objekt
        return arc.intersect_circle(self, pickdist, include_end)
        
    def split_OLD(self, x, pickdist):
        """split this circle at intersection point(s) with line or at points"""
        if isinstance(x, Line):
            line = x
            points = sorted(self.intersect(line, pickdist))
        else:
            points = sorted([p for p in x
                             if (np.isclose(np.linalg.norm(p), self.radius) and
                                 self.is_inside(np.arctan2(p[1], p[0]),
                                                pickdist))])
            
        elements = []
        if len(points) > 1:
            z1 = np.asarray(points[0]) - np.asarray(self.center)
            z2 = np.asarray(points[1]) - np.asarray(self.center)
            startangle = np.arctan2(z1[1], z1[0])
            endangle = np.arctan2(z2[1], z2[0])
            elements.append(Arc(
                Element(center=self.center,
                        radius=self.radius,
                        start_angle=startangle*180/np.pi,
                        end_angle=endangle*180/np.pi)))
        return elements

    def XXX_split(self, e):
        if isinstance(e ,Line):
            return self.split_line(e)
            
        return e
        
    def __str__(self):
        return "Circle c={}, r={}".format(self.center, self.radius)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self == other

    def __hash__(self):
        """Override the default hash behavior (that returns the id or the object)"""
        return hash(tuple(sorted(self.__dict__.items())))

#############################
#        Arc (Shape)        #
#############################

class Arc(Circle):
    """a counter clockwise segment of a circle with start and end point"""
    def __init__(self, e):
        super(self.__class__, self).__init__(e)
        self.startangle = e.start_angle/180*np.pi
        self.endangle = e.end_angle/180*np.pi
        if self.endangle < self.startangle:
            if self.endangle < 0:
                self.endangle += 2*np.pi
            elif self.startangle < 0:
                self.startangle += 2*np.pi
            else:
                self.endangle -= 2*np.pi

        logger.debug("%f -- %f (C %f, %f R %f)",
                     self.startangle, self.endangle,
                     self.center[0], self.center[1], self.radius)
        self.p1 = (self.center[0] + e.radius*np.cos(self.startangle),
                   self.center[1] + e.radius*np.sin(self.startangle))
        self.p2 = (self.center[0] + e.radius*np.cos(self.endangle),
                   self.center[1] + e.radius*np.sin(self.endangle))
                            
    def render(self, renderer, color='blue', with_nodes=False):
        renderer.arc(self.startangle, self.endangle,
                     self.center, self.radius, color)
        if with_nodes:
            renderer.point(self.p1, 'ro', color)
            renderer.point(self.p2, 'ro', color)

    def center_of_connection(self):
        s = self.startangle
        d = self.endangle - s
        if d > 2*np.pi:
            d -= 2*np.pi
        x, y = self(s + d/2)
        return (round(x, ndec), round(y, ndec))

    def length(self):
        """returns length of this arc"""
        d = abs(self.endangle - self.startangle)
        if d > 2*np.pi:
            d -= 2*np.pi
        return self.radius*abs(d)

    def range(self, step=1.0):
        """returns evenly spaced values"""
        num = max(self.length()/step, 1)
        astep = self.length()/num/self.radius
        s = self.startangle
        d = self.endangle - s
        if d > 2*np.pi:
            d -= 2*np.pi
        alpha = np.arange(s, s+d, astep)
        return self(alpha)
    
    def __call__(self, alpha):
        """returns x,y coordinates of angle"""
        return (self.center[0] + self.radius*np.cos(alpha),
                self.center[1] + self.radius*np.sin(alpha))

    def intersect(self, line, pickdist=1e-3, include_end=False):
        """calculate intersection of this arc with line"""
        logger.debug(self)
        points = super(Arc, self).intersect(line, pickdist)
        ip = []
        logger.debug("points %s", points)
        for p in points:
            z = np.asarray(p) - np.asarray(self.center)
            alpha = np.arctan2(z[1], z[0])
            if self.is_inside(alpha, pickdist, include_end):
                ip.append(p)
        return ip
    
    def split_OLD(self, x, pickdist):
        """split this arc at intersection point(s) with line"""
        if isinstance(x, Line):
            line = x
            points = sorted(self.intersect(line, pickdist))
        else:
            points = sorted([p for p in x
                             if (np.isclose(np.linalg.norm(p), self.radius) and
                                 self.is_inside(np.arctan2(p[1], p[0]),
                                                pickdist))])

        elements = []
        for p in points:
            z = np.asarray(p) - np.asarray(self.center)
            alpha = np.arctan2(z[1], z[0])
            if self.is_inside(alpha, pickdist, include_end=True):
                elements.append(Arc(
                    Element(center=self.center,
                            radius=self.radius,
                            start_angle=self.startangle*180/np.pi,
                            end_angle=alpha*180/np.pi)))
                elements.append(Arc(
                    Element(center=self.center,
                            radius=self.radius,
                            start_angle=alpha*180/np.pi,
                            end_angle=self.endangle*180/np.pi)))
        return elements

    def intersect_line(self, line, pickdist=1e-3, include_end=False):
        """ Von einem Arc-Objekt und einem Line-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben
        """
        points = super(Arc, self).intersect_line(line, pickdist, include_end)
        
        # Die Funktion der Basis Circle hat die möglichen Punkte bestimmt.
        # Nun wird geprüft, ob sie auf dem Kreissegment liegen.
        remaining_points = []
        for p in points:
            if self.is_point_inside(p, pickdist, include_end):
                remaining_points.append(p)
        return remaining_points

    def intersect_arc(self, arc, pickdist=1e-3, include_end=False):
        """ Von zwei Arc-Objekten werden die Schnittpunkte bestimmt und in
            einer Liste ausgegeben.
        """
        assert(isinstance(arc, Arc))
        
        points = self.intersect_circle(arc, pickdist, include_end)

        # Nun wird geprüft, ob die Punkte auch auf dem arc-Kreissegment liegen.
        # Dieses wurde bis jetzt als Kreis betrachtet.
        remaining_points = []
        for p in points:
            if arc.is_point_inside(p, pickdist, include_end):
                remaining_points.append(p)
        return remaining_points

    def intersect_circle(self, circle, pickdist=1e-3, include_end=False):
        """ Von einem Arc-Objekt und einem Circle-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben
        """
        if points_are_close(self.center, circle.center, pickdist):
            if np.isclose(self.radius, circle.radius):
                if include_end:
                    return [self.p1, self.p2]
            # Wenn bei gleichem Mittelpunkt der Radius abweicht, gibt es sicher
            # keine Schnittpunkt
            return []
            
        points = super(Arc, self).intersect_circle(circle, pickdist, include_end)
        
        # Die Schnittpunkte von zwei Circle-Objekten sind bestimmt. Nun wird
        # geprüft, ob sie auf dem Kreissegment liegen.
        remaining_points = []
        for p in points:
            if self.is_point_inside(p, pickdist, include_end):
                remaining_points.append(p)
        return remaining_points

    def is_point_inside(self, p, pickdist, include_end=False):
        """ Die Funktion prüft, ob der Punkt p auf dem Kreissegment liegt.
        """
        if points_are_close(p, self.p1, pickdist, 1e-4):
            return include_end
        elif points_are_close(p, self.p2, pickdist, 1e-4):
            return include_end
        elif points_are_close(self.p1, self.p2, pickdist, 1e-4):
            return False

        alpha_p1 = alpha_line(self.center, self.p1)
        alpha_p2 = alpha_line(self.center, self.p2)
        alpha_p = alpha_line(self.center, p)
        alpha_inside = is_inside(alpha_p1, alpha_p2, alpha_p)

        dist = distance(self.p1, self.p2)
        dist_p1 = distance(self.p1, p)
        dist_p2 = distance(self.p2, p)
        dist_inside = not (dist_p1 > dist or dist_p2 > dist)

        if dist_inside == False and alpha_inside == True:
            print("FATAL BUG")
            print("   p1={}".format(self.p1))
            print("   pm={}".format(p))
            print("   p2={}".format(self.p2))
            print("   inside dist={}, alpha={}".format(dist_inside, alpha_inside))
            assert(False)
            
        return alpha_inside
        
    def is_inside(self, alpha, pickdist, include_end=False):
        """ returns True if alpha is between start and end angle
        """
        return is_inside(self.startangle, self.endangle, alpha)

    def transform(self, T, **kwargs):
        super(Arc, self).transform(T)
        p1, p2 = ((self.p1[0]-self.center[0],
                   self.p1[1]-self.center[1]),
                  (self.p2[0]-self.center[0],
                   self.p2[1]-self.center[1]))

        if kwargs.get('reflect', False):
            self.p1, self.p2 = self.p2, self.p1
            self.endangle = np.arctan2(p1[1], p1[0])
            self.startangle = np.arctan2(p2[1], p2[0])
        else:
            self.startangle = np.arctan2(p1[1], p1[0])
            self.endangle = np.arctan2(p2[1], p2[0])
        return self

    def minmax(self):
        """ Die Funktion bestimmt das Minimum und Maximum auf der x- und der
            y-Achse (return [<min-x>, <max-x>, <min-y>, <max-y>])
        """
        mm = [min(self.p1[0], self.p2[0]), max(self.p1[0], self.p2[0]),
              min(self.p1[1], self.p2[1]), max(self.p1[1], self.p2[1])]

        p = [self.center[0]-self.radius, self.center[1]]
        if p[0] < mm[0]:
            a = alpha_line(self.center, p)
            if self.is_inside(a, 0.00001):
                mm[0] = p[0]
            
        p = [self.center[0]+self.radius, self.center[1]]            
        if p[0] > mm[1]:
            a = alpha_line(self.center, p)
            if self.is_inside(a, 0.00001):
                mm[1] = p[0]

        p = [self.center[0], self.center[1]-self.radius]
        if p[1] < mm[2]:
            a = alpha_line(self.center, p)
            if self.is_inside(a, 0.00001):
                mm[2] = p[1]
            
        p = [self.center[0], self.center[1]+self.radius]
        if p[1] > mm[3]:
            a = alpha_line(self.center, p)
            if self.is_inside(a, 0.00001):
                mm[3] = p[1]
        return mm

    def minmax_from_center(self, center):
        """ Die Funktion ermittelt den minimalen und maximalen Abstand vom Center 
        """
        d = distance(center, self.center)
        if np.isclose(d, 0.0):
            return (self.radius, self.radius)
            
        angle = alpha_line(center, self.center)
        dist_min = abs(d - self.radius)
        dist_max = d + self.radius
        
        pmax = point(center, d + self.radius, angle)
        alpha_pmax = alpha_line(self.center, pmax)
        if not self.is_inside(alpha_pmax, 1e-08):
            dist_max = max(distance(center, self.p1), distance(center, self.p2))

        pmin = point(center, d - self.radius, angle)
        alpha_pmin = alpha_line(self.center, pmin)

        if not self.is_inside(alpha_pmin, 1e-08):
            dist_min = min(distance(center, self.p1), distance(center, self.p2))

        return (dist_min, dist_max)

    def get_nodes(self):
        """ Die Funktion liefert eine Liste von virtuellen Nodes, welche man
            zum Rechnen der convex_hull() benötigt.
        """
        return (p for p in points_on_arc(self.center, self.radius,
                                         self.startangle, self.endangle))
        
    def __str__(self):
        return "Arc c={}, r={} start={}, end={}, p1={}, p2={}".format(self.center,
                                                                self.radius, self.startangle,
                                                                self.endangle,
                                                                self.p1, self.p2)
                            
    def __eq__(self, other):
        """Override the default Equals behavior"""
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self == other

    def __hash__(self):
        """Override the default hash behavior (that returns the id or the object)"""
        return hash(tuple(sorted(self.__dict__.items())))

#############################
#        Line (Shape)       #
#############################

class Line(Shape):
    """straight connection between start and end point"""
    def __init__(self, e):
        self.p1 = e.start[0], e.start[1]
        self.p2 = e.end[0], e.end[1]

    def render(self, renderer, color='blue', with_nodes=False):
        renderer.line(self.p1, self.p2, color)
        if with_nodes:
            renderer.point(self.p1, 'ro', color)
            renderer.point(self.p2, 'ro', color)

    def center_of_connection(self):
        if np.isclose(self.dx(), 0):
            return (round(self.p1[0], ndec),
                    round(self.dy()/2. + self.ymin(), ndec))
        x = self.p1[0] + self.dx()/2.
        return (round(x, ndec), round(self(x), ndec))

    def length(self):
        return np.sqrt(self.dx()**2 + self.dy()**2)
    
    def range(self, step=1.0):
        """returns evenly spaced values"""
        num = max(int(self.length()/step), 1)
        if np.isclose(self.dx(), 0):
            if np.isclose(self.dy(), 0):
                return ((self.xmin(),), (self.ymin(),))
            y = np.arange(self.ymin(), self.ymax(), self.length()/num)
            return np.array((self.xmin()*np.ones(len(y)), y))
        xstep = np.sqrt((self.length()/num)**2/(1 + self.dy()**2/self.dx()**2))
        x = np.arange(self.xmin(), self.xmax(), xstep)
        return (x, self(x))
    
    def __call__(self, x):
        """returns y coordinate of x"""
        if np.isclose(self.dx(), 0):
            return float('nan')
        return self.dy()/self.dx()*(x - self.p1[0]) + self.p1[1]

    def intersect(self, line, pickdist=1e-3, include_end=False):
        """returns intersection point with line"""
        point = []
        logger.debug(self)
        if np.isclose(self.dx(), 0):
            if np.isclose(line.dx(), 0):
                return []
            point = (self.p1[0], line(self.p1[0]))
        elif np.isclose(line.dx(), 0.0):
            point = (line.p1[0], self(line.p1[0]))
        else:
            mself, mline = (self.dy()/self.dx(), line.dy()/line.dx())
            dm = mself - mline
            dc = (line.p1[1] - mline*line.p1[0]) - \
                 (self.p1[1] - mself*self.p1[0])
            if np.isclose(dm, 0.0):
                # check start and end point
                if la.norm(np.asarray(line.p1) - self.p1) <= pickdist or \
                   la.norm(np.asarray(line.p2) - self.p1) <= pickdist:
                    return (self.p1,)
                if la.norm(np.asarray(line.p1) - self.p2) <= pickdist or \
                   la.norm(np.asarray(line.p2) - self.p2) <= pickdist:
                    return (self.p2,)
                return ()
            x = dc/dm
            point = (x, self(x))
        if self.is_inside(point, pickdist, include_end):
            return (point,)
        return []
        
    def intersect_line(self, line, pickdist=1e-3, include_end=False):
        """ Von zwei Line-Objekten wird der Schnittpunkt bestimmt und in
            einer Liste ausgegeben.
        """
        point = []
        m_L1 = self.m()
        m_L2 = line.m()
        if m_L1 == None:
            if m_L2 == None:
                return []
            else:
                y = line_n([line.p1[0]-self.p1[0], line.p1[1]], m_L2)
                point = [self.p1[0], y]
        else:
            if m_L2 == None:
                y = line_n([self.p1[0]-line.p1[0], self.p1[1]], m_L1)
                point = [line.p1[0], y]
            else:
                if np.isclose(m_L1, m_L2):
                    return []
                else:
                    point = lines_intersect_point(m_L1, self.n(m_L1), m_L2, line.n(m_L2))
                    
        if line.is_inside(point, pickdist, include_end):
            if self.is_inside(point, pickdist, include_end):
                return [point]
        return []

    def intersect_arc(self, arc, pickdist=1e-3, include_end=False):
        """ Von einem Line-Objekt und einem Arc-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben.
        """
        return arc.intersect_line(self, pickdist, include_end)

    def intersect_circle(self, circle, pickdist=1e-3, include_end=False):
        """ Von einem Line-Objekt und einem Circle-Objekt werden die
            Schnittpunkte bestimmt und in einer Liste ausgegeben.
        """
        return circle.intersect_line(self, pickdist, include_end)
        
    def split_OLD(self, line, pickdist):
        """split this line at intersection point with line"""
        point = self.intersect(line, pickdist)
        if point and self.is_inside(point[0], pickdist):
            logger.info("Line intersection %s", point)
            return (Line(Element(start=self.p1, end=point[0])),
                    Line(Element(start=point[0], end=self.p2)))
        return []
    
    def is_inside(self, point, pickdist, include_end=False):
        """ returns True if point is between start and end point
        """
        logger.debug("Arc::is_inside: %s in (%s, %s)", point, self.p1, self.p2)
        
        if points_are_close(point, self.p1, pickdist):
            return include_end
        if points_are_close(point, self.p2, pickdist):
            return include_end
            
        m1 = line_m(self.p1, point)
        m2 = line_m(point, self.p2)
        if m1 == None or m2 == None:
            if m1 != None or m2 != None:
                return False
        elif not np.isclose(m1, m2, pickdist):
            return False
        
        length = distance(self.p1, self.p2)
        dist_p1 = distance(self.p1, point)
        dist_p2 = distance(self.p2, point)
        if dist_p1 > length or dist_p2 > length:
            return False
        return True
         
    def minmax(self):
        """ Die Funktion bestimmt das Minimum und Maximum auf der x- und der
            y-Achse (return [<min-x>, <max-x>, <min-y>, <max-y>])
        """
        return [min(self.p1[0], self.p2[0]), max(self.p1[0], self.p2[0]),
                min(self.p1[1], self.p2[1]), max(self.p1[1], self.p2[1])]
                
    def minmax_from_center(self, center):
        """ Die Funktion ermittelt den minimalen und maximalen Abstand vom Center 
        """
        dist_max = max(distance(center, self.p1), distance(center, self.p2))
        dist_min = min(distance(center, self.p1), distance(center, self.p2))
        return (dist_min, dist_max)
        
    def get_nodes(self):
        """ Die Funktion liefert eine Liste von virtuellen Nodes, welche man
            zum Rechnen der convex_hull() benötigt.
        """
        return (self.p1, self.p2)

    def __str__(self):
        return "Line p1={}, p2={}".format(self.p1, self.p2)

    def __eq__(self, other):
        """Override the default Equals behavior"""
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self == other

    def __hash__(self):
        """Override the default hash behavior (that returns the id or the object)"""
        return hash(tuple(sorted(self.__dict__.items())))

#############################
#       Point (Shape)       #
#############################
    
class Point(Shape):
    """ Das Point-Objekt wird nur für die Ausgabe im Plot verwendet """
    def __init__(self, p):
        self.p1 = p

    def render(self, renderer):
        renderer.point(self.p1)


def is_connected(a, b):
    if a == b:
        return False
    return a[0] == b[0] or a[1] == b[0] or a[0] == b[1] or a[1] == b[1]


def single_path(edges):
    """sort edges to create a sequence of connected edges"""
    if edges:
        p = [edges[0]]
        remaining = edges[1:]
        while remaining:
            edges = remaining[:]
            remaining = []
            for e in edges:
                if is_connected(p[-1], e):
                    p.append(e)
                elif is_connected(p[0], e):
                    p.insert(0, e)
                else:
                    remaining.append(e)
        return p
    return []

#############################
#           Motor           #
#############################

class Motor(object):
    def __init__(self, center, radius, startangle=0, endangle=0):
        self.center = center
        self.radius = radius
        self.startangle = startangle
        self.endangle = endangle
        self.part = self.part_of_circle()
        self.airgaps = []
        self.airgap_radius = 0.0
        
        print("Motor: Center=({},{}), Radius={}, Start={}, End={}"
            .format(self.center[0], self.center[1], self.radius, self.startangle, self.endangle))

    def is_a_motor(self):
        return self.radius > 0.0
        
    def is_in_middle(self):
        return self.radius > 0.0 and points_are_close(self.center, [0.0, 0.0], 1e-8)

    def is_full(self):
        return self.radius > 0.0 and self.startangle == 0.0 and self.endangle == 0.0
        
    def is_half_up(self):
        return self.radius > 0.0 and \
               np.isclose(self.startangle, 0.0, 1e-8) and \
               (np.isclose(self.endangle, np.pi, 1e-8) or \
                np.isclose(self.endangle, -np.pi, 1e-8))

    def is_half_down(self):
        return self.radius > 0.0 and \
               (np.isclose(self.endangle, np.pi, 1e-8) or \
                np.isclose(self.endangle, -np.pi, 1e-8)) and\
               np.isclose(self.startangle, 0.0, 1e-8)

    def is_half_left(self):
        return self.radius > 0.0 and \
               np.isclose(self.startangle, np.pi/2, 1e-8) and \
               np.isclose(self.endangle, -np.pi/2, 1e-8)

    def is_half_right(self):
        return self.radius > 0.0 and \
               np.isclose(self.startangle, -np.pi/2, 1e-8) and \
               np.isclose(self.endangle, np.pi/2, 1e-8)
               
    def is_half(self):
        return self.is_half_up() or self.is_half_down() or \
               self.is_half_left() or self.is_half_right()

    def is_quarter(self):
        return self.radius > 0.0 and \
               np.isclose(alpha_angle(self.startangle, self.endangle), np.pi/2)
        
    def move_to_middle(self, geom):
        if not self.is_in_middle():
            if self.radius > 0.0:
                offset = [-(self.center[0]), -(self.center[1])]
                geom.move(offset)
                self.set_center(geom, 0.0, 0.0)
                geom.clear_schnitt()
                return True
        return False

    def set_center(self, geom, x, y):
        self.center[0] = x
        self.center[1] = y
        geom.center = [x, y]

    def set_radius(self, radius):
        self.radius = radius
        
    def copy(self, geom, startangle, endangle, airgap=False, inside=True):
        if airgap and self.airgap_radius > 0.0:
            if inside:
                clone = geom.copy_shape( self.center, self.radius, startangle, endangle, 0.0, self.airgap_radius)
            else:
                clone = geom.copy_shape( self.center, self.radius, startangle, endangle, self.airgap_radius, self.radius)
        else:
            clone = geom.copy_shape( self.center, self.radius, startangle, endangle, 0.0, self.radius+50)

        if not np.isclose(normalise_angle(startangle), normalise_angle(endangle), 0.0):
            start_p = point(self.center, self.radius+5, startangle)
            start_line = Line(Element(start=self.center, end=start_p))
            clone.add_schnitt(start_line)
            
            end_p = point(self.center, self.radius+5, endangle)
            end_line = Line(Element(start=self.center, end=end_p))
            clone.add_schnitt(end_line)
        
        if self.radius > 0.0:
            circ = Circle(Element(center=self.center, radius=self.radius+1))
            clone.add_schnitt(circ)
        return clone

    def rotate_to(self, geom, new_startangle):
        if points_are_close(self.center, [0.0, 0.0]):
            angle = new_startangle - self.startangle
            geom.rotate(angle)
            self.startangle -= angle
            self.endangle -= angle
        else:
            print("rotate nur mit Center(0,0) möglich")
        
    def airgap(self, geom):
        if self.radius > 0.0:
            self.airgaps = geom.detect_airgaps(self.center, self.startangle, self.endangle)
            if len(self.airgaps) > 0:
                num_airgaps = 0
                for g in self.airgaps:
                    print("Airgap zwischen {}".format(g))
                    gap_radius = round((g[0]+g[1])/2.0, 6)
                    circle = Circle(Element(center=self.center, radius=gap_radius))
                    if geom.is_airgap(self.center, self.radius, self.startangle, self.endangle, circle):
                        print("Airgap with radius {}".format(gap_radius))
                        geom.airgaps.append(circle)
                        num_airgaps += 1
                        self.airgap_radius = gap_radius
                    else:
                        geom.airgaps.append(circle)
                        print("DESASTER: No Airgap with radius {}".format(gap_radius))

                if num_airgaps == 1:
                    print("One Airgap found !!! {}".format(self.airgap_radius))
                else:
                    self.airgap_radius = 0.0
                    
    def has_airgap(self):
        return self.airgap_radius > 0

    def part_of_circle(self):
        return part_of_circle(self.startangle, self.endangle)

    def repair_hull(self, geom):
        if self.part > 1:
            geom.repair_hull_line(self.center, self.startangle)
            geom.repair_hull_line(self.center, self.endangle)
        
    def find_symmetrie(self, geom):
        if self.radius > 0.0:
            geom.find_symmetrie(self.center, self.radius)

#############################
#         Geometrie         #
#############################

class Geometry(object):
    """collection of connected shapes"""
    def __init__(self, elements=[], pickdist=1e-3, adjust=True):
        self._name = ''
        self.mirror_axis = None
        self.corners = []
        self.center = []
        self.schnitt = []
        self.symlines = []
        self.airgaps = []
        self.g = nx.Graph()
        self.pickdist = pickdist
        i = 0
#        for e in add_and_split(elements, self.pickdist):
        for e in elements:
            if e:
                e.id = i
                n = self.find_nodes(e.start(), e.end())
                logger.debug("%d %s", i, n)
                try:
                    add_or_join(self.g, n[0], n[1], e, pickdist)
                except Exception as ex:
                    logger.warn("EXCEPTION %s", ex)
                    if e:  # must be a circle
                        self.g.add_node(e.center, object=e)
            i += 1

        if self.g.nodes():
            offset = np.min(np.array(self.g.nodes()), axis=0)
            logger.info("xmin, ymin %s", offset)
            if not np.isclose(offset, (0, 0)).all() and adjust:
                self.move(-offset)
            # set corners
            self.corners = sorted(find_corners(self.g.nodes()))
            logger.info("Corners %s", self.corners)
            # set section angle
            self.alpha = np.pi/round(np.pi/self.argmax())

            self.set_diameters()

            shaft = self.shaft()
            if shaft:
                self.diameters = (round(shaft, 3),
                                  self.diameters[-1])

    def shaft(self):
        """returns shaft diameter if any"""
        radius = []
        for a in self.arcs():
            if np.isclose(a.center, 0.0).all():
                logger.info("Shaft candidate (%f, %f) %r",
                            a.startangle/np.pi*180,
                            a.endangle/np.pi*180,
                            a.radius)
                if 0 < a.radius < self.diameters[-1]/2:
                    radius.append(a.radius)
        return 2*sorted(radius)[0] if radius else 0
    
    def set_diameters(self):
        try:
            min, max = self.min(), self.max()
            logger.info("Min %s max %s", min, max)
            self.diameters = (round(2*np.sqrt(min[0]**2 + min[1]**2), 3),
                              round(2*np.sqrt(max[0]**2 + max[1]**2), 3))
        except:
            self.diameters = []
            
    def move(self, offset):
        """moves all objects"""
        for e in self.g.edges(data=True):
            e[2]['object'].move(offset)
        for c in self.circles():
            c.move(-offset)
        mapping = {n: (round(n[0]+offset[0], 3), round(n[1]+offset[1], 3))
                   for n in self.g.nodes()}
        nx.relabel_nodes(self.g, mapping, copy=False)
        
    def rotate(self, alpha):
        """rotates all objects by angle alpha"""
        T = np.array(((np.cos(alpha), -np.sin(alpha)),
                      (np.sin(alpha), np.cos(alpha))))
        for e in self.g.edges(data=True):
            e[2]['object'].transform(T)
        for c in self.circles():
            c.transform(T)
        rotnodes = np.dot(T, np.asarray(self.g.nodes()).T).T.tolist()
        mapping = {n: (round(r[0], 3), round(r[1], 3))
                   for n, r in zip(self.g.nodes(), rotnodes)}
        nx.relabel_nodes(self.g, mapping, copy=False)
        
    def scale(self, factor):
        """scales all objects"""
        for e in self.g.edges(data=True):
            e[2]['object'].scale(factor)
        for c in self.circles():
            c.scale(factor)
        mapping = {n: (round(factor*n[0], 3), round(factor*n[1], 3))
                   for n in self.g.nodes()}
        nx.relabel_nodes(self.g, mapping, copy=False)
        self.diameters = tuple([factor*d for d in self.diameters])

    def find_nodes0(self, *points):
        """return closest nodes to points in arg within pickdist"""
        return [tuple([round(int(x/self.pickdist+0.5)*self.pickdist, ndec)
                       for x in p])
                for p in points]
        
    def find_nodes(self, *points, **kwargs):
        """return closest nodes to points in arg within pickdist"""
        n = []
        nodes = kwargs.get('g', self.g).nodes()
        if nodes:
            anodes = np.asarray(nodes)
            for p in points:
                # la.norm on numpy below 1.8 does not accept axis
                c = anodes - p
                dist = np.sqrt(np.einsum('ij, ij->i', c, c))
                # dist = la.norm(np.asarray(nodes) - p, axis=1)
                idx = dist.argmin()
                if dist[idx] < self.pickdist:
                    n.append(nodes[idx])
                else:
                    n.append((round(p[0], ndec), round(p[1], ndec)))
        else:
            return [(round(p[0], ndec), round(p[1], ndec)) for p in points]
        return n
    
    def polyline_paths(self):
        """return lists of line paths"""
        paths = []
        g = self.g.copy()
        g.remove_edges_from(self.arcs())
        while g.number_of_edges() > 0:
            p = single_path(g.edges())
            g.remove_edges_from(p)
            # rearrange sequence to make it contiguous:
            px = [p[0]]
            for e in p[1:]:
                if px[-1][1] == e[1]:
                    px.append(tuple(reversed(e[:2])))
                elif px[0][0] == e[0]:
                    px.insert(0, tuple(reversed(e[:2])))
                elif px[0][0] == e[1]:
                    px.insert(0, e[:2])
                else:
                    px.append(e[:2])
            # only keep end point of each edge
            px[1:] = [x[-1] for x in px]
            px[0] = px[0][0]
            paths.append(px)
        return paths

    def number_of_nodes(self):
        """returns the number of nodes in graph"""
        return self.g.number_of_nodes()
    
    def number_of_edges(self):
        """return the number of edges in graph"""
        return self.g.number_of_edges()
        
    def elements(self, type):
        """return lists of objects"""
        return [e[2]['object'] for e in self.g.edges(data=True)
                if isinstance(e[2]['object'], type)]
        
    def arcs(self):
        """return lists of arcs"""
        return self.elements(Arc)
    
    def lines(self):
        """return lists of lines"""
        return self.elements(Line)

    def circles(self):
        """return list of circle nodes"""
        return [n[1]['object'] for n in self.g.nodes(data=True)
                if n[1] and isinstance(n[1]['object'], Circle)]

    def virtual_nodes(self):
        nodes = []
        for e in self.elements(Shape):
            nodes += e.get_nodes()
        return nodes

    def angle_nodes(self, center, angle, rtol, atol):
        nodes = []
        for n in self.g.nodes():
            if points_are_close(center, n, rtol, atol):
                # Da gibt es leinen brauchbareb Winkel
                nodes.append(n)
            elif np.isclose(angle, alpha_line(center, n), rtol, atol):
                nodes.append(n)
        return nodes

    def repair_hull_line(self, center, angle):
        # Bei der Suche sind wir bei der pickdist für isclose() tolerant,
        # sonst finden wir einige Nodes nicht
        rtol = 1e-2
        atol = 1e-4
        
        corners = [Corner(center,c) for c in self.angle_nodes(center, angle, rtol, atol)]
        if len(corners) < 2:
            # Ohne genügend Corners ist die Arbeit sinnlos
            return
        corners.sort()

        for p1, p2 in self.g.edges():
            for c in corners:
                if c.is_equal(p1):
                    if not (points_are_close(center, p2, rtol, atol) or \
                            np.isclose(angle, alpha_line(center, p2), rtol, atol)):
                        c.set_keep_node()
                elif c.is_equal(p2):
                    if points_are_close(center, p1, rtol, atol) or \
                       np.isclose(angle, alpha_line(center, p1), rtol, atol):
                        self.g.remove_edge(p1, p2)
                    else:
                        c.set_keep_node()

        for c in corners:
            if not c.keep_node():
                self.g.remove_node(c.point())

        # Weil uns unnötige Corners abhanden gekommen sind, bilden wir die
        # Liste neu.
        corners = [Corner(center, c) for c in self.angle_nodes(center, angle, rtol, atol)]
        if len(corners) > 1:
            corners.sort()
            p1 = corners[0].point()
            for c in corners[1:]:
                p2 = c.point()
                self.g.add_edge(p1, p2, object=Line(Element(start=p1, end=p2)))
                p1 = p2

    def remove_corners(self, g):
        """removes the corner nodes with their edges"""
        corners = [n for n in self.find_corners(g.nodes()) if g.degree(n) < 3]
        logger.debug("removing corners %s", corners)
        g.remove_nodes_from(corners)
           
    def area_nodes(self, incl_bnd=False):
        """return nodes that enclose areas (eg. slots)"""
        g = self.g.copy()
        if not incl_bnd:
            self.remove_corners(g)
        return nx.cycle_basis(g)
            
    def areas(self, incl_bnd=False):
        """return sorted paths that enclose areas (eg. slots)"""
        areas = []
        for nodes in self.area_nodes(incl_bnd):
            edges = [e for e in self.g.edges(data=True)
                     if e[0] in nodes and e[1] in nodes]
            areas.append([e[2]['object'] for e in single_path(edges)])
        return areas
        
    def remove_areas(self, incl_bnd=False):
        """remove all edges that are connected to area nodes"""
        for nodes in self.area_nodes(incl_bnd):
            edges = [e for e in self.g.edges()
                     if e[0] in nodes and e[1] in nodes]
            self.g.remove_edges_from(edges)

    def find_corners(self, nodes):
        """find corners of nodes"""
        if nodes:
            minX, minY = np.amin(nodes, axis=0)
            maxX, maxY = np.amax(nodes, axis=0)
            print("find_corners: minX={}, minY={}, maxX={}, maxY={}".format(minX, minY, maxX, maxY))
            return [n for n in nodes
                    if n[0] in (minX, maxX) or n[1] in (minY, maxY)]
        return []
     
    def find_nearest_node(self, ref, incl_bnd=True):
        """return nearest node to ref node of all areas"""
        n = []
        for nodes in self.area_nodes(incl_bnd):
            n += nodes
        indx = np.argmin(np.array([np.linalg.norm(z)
                                   for z in np.array(n) - np.asarray(ref)]))
        return n[indx]

    def copy_line(self, center, radius, start_angle, end_angle,
                  start_line, end_line, inner_circle, outer_circle, e):
        """ Die Funktion kopiert die Teile einer Linie, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Line))
        if is_same_angle(start_angle, end_angle):
            points = inner_circle.intersect_line(e, self.pickdist, False) + \
                     outer_circle.intersect_line(e, self.pickdist, False) + \
                     [e.p2]
        else:
            points = e.intersect_line(start_line, self.pickdist, False) + \
                     e.intersect_line(end_line, self.pickdist, False) + \
                     inner_circle.intersect_line(e, self.pickdist, False) + \
                     outer_circle.intersect_line(e, self.pickdist, False) + \
                     [e.p2]

        new_elements = []
        sorted_points = []

        for p in points:
            dist = distance(e.p1, p)
            sorted_points.append((dist ,p))
        sorted_points.sort()

        p1 = e.p1
        for x,p2 in sorted_points:
            pm = middle_of_line(p1, p2)
            if is_inside_region(pm, center, inner_circle.radius, outer_circle.radius,
                                start_angle, end_angle):
                new_elements.append(Line(Element(start=p1, end=p2)))
            p1 = p2
            
        return new_elements
        
    def copy_arc(self, center, radius, start_angle, end_angle,
                 start_line, end_line, inner_circle, outer_circle, e):

        """ Die Funktion kopiert die Teile eines Kreissegments, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Arc))
        if is_same_angle(start_angle, end_angle):
            points = inner_circle.intersect_arc(e, self.pickdist, False) + \
                     outer_circle.intersect_arc(e, self.pickdist, False) + \
                     [e.p2]
        else:
            points = e.intersect_line(start_line, self.pickdist, False) + \
                     e.intersect_line(end_line, self.pickdist, False) + \
                     inner_circle.intersect_arc(e, self.pickdist, False) + \
                     outer_circle.intersect_arc(e, self.pickdist, False) + \
                     [e.p2]
           
        new_elements = []
        sorted_points = []
        
        alpha_start = alpha_line(e.center, e.p1)
        
        for p in points:
            alpha_next = alpha_line(e.center, p)
            if less_equal(alpha_next, alpha_start):
                alpha_next += 2*np.pi
            sorted_points.append((alpha_next ,p))
            alpha_start = alpha_next
        sorted_points.sort()

        p1 = e.p1
        alpha_start = alpha_line(e.center, e.p1)
        for x,p2 in sorted_points:
            alpha_end = alpha_line(e.center, p2)
            pm = middle_of_arc(e.center, e.radius, p1, p2)
            if is_inside_region(pm, center, inner_circle.radius, outer_circle.radius,
                                start_angle, end_angle):
                new_elements.append(Arc(Element(center=e.center, radius=e.radius,
                                                start_angle=alpha_start*180/np.pi,
                                                end_angle=alpha_end*180/np.pi)))
            
            alpha_start = alpha_end
            p1 = p2
        return new_elements

    def copy_circle(self, center, radius, start_angle, end_angle,
                    start_line, end_line, inner_circle, outer_circle, e):
        """ Die Funktion kopiert die Teile eines Kreises, welche sich in der
            durch die Parameter definierten Teilkreisfläche befinden.
        """
        assert(isinstance(e, Circle))
        if is_same_angle(start_angle, end_angle):
            points = inner_circle.intersect_circle(e, self.pickdist, False) + \
                     outer_circle.intersect_circle(e, self.pickdist, False)
        else:
            points = e.intersect_line(start_line, self.pickdist) + \
                     e.intersect_line(end_line, self.pickdist) + \
                     inner_circle.intersect_circle(e, self.pickdist, False) + \
                     outer_circle.intersect_circle(e, self.pickdist, False)
                 
        new_elements = []
        if len(points) < 2:
            if is_inside_region(e.center, center, inner_circle.radius, outer_circle.radius,
                                start_angle, end_angle):
               new_elements.append(Circle(Element(center=e.center, radius=e.radius)))
            return new_elements
        
        sorted_points = []
        for p in points:
            alpha_p = alpha_line(e.center, p)
            sorted_points.append((alpha_p, p))
        sorted_points.sort()

        x,px = sorted_points[0]
        del sorted_points[0]
        p1 = px
        alpha_start = alpha_line(e.center, p1)
        for x,p2 in sorted_points:
            alpha_end = alpha_line(e.center, p2)
            pm = middle_of_arc(e.center, e.radius, p1, p2)
            if is_inside_region(pm, center, inner_circle.radius, outer_circle.radius,
                                start_angle, end_angle):
                new_elements.append(Arc(Element(center=e.center, radius=e.radius,
                                                start_angle=alpha_start*180/np.pi,
                                                end_angle=alpha_end*180/np.pi)))
            
            alpha_start = alpha_end
            p1 = p2

        alpha_end = alpha_line(e.center, px)
        pm = middle_of_arc(e.center, e.radius, p1, px)
        if is_inside_region(pm, center, inner_circle.radius, outer_circle.radius,
                            start_angle, end_angle):
            new_elements.append(Arc(Element(center=e.center, radius=e.radius,
                                            start_angle=alpha_start*180/np.pi,
                                            end_angle=alpha_end*180/np.pi)))
        return new_elements
      
    def copy_shape(self, center, radius, startangle, endangle, inner_radius, outer_radius):
        """ Die Funktion kopiert die Teile von Shape-Objekten, welche sich in
            der durch die Parameter definierten Teilkreisfläche befinden.
        """
        if is_same_angle(startangle, endangle):
            start_line = Line(Element(start=center, end=point(center, radius+1, startangle)))
            end_line = Line(Element(start=center, end=point(center, radius+1, startangle)))
        else:
            start_line = Line(Element(start=center, end=point(center, radius+1, startangle)))
            end_line = Line(Element(start=center, end=point(center, radius+1, endangle)))

        if np.isclose(normalise_angle(startangle), normalise_angle(endangle), 0.0):
            inner_circle = Circle(Element(center=center, radius=inner_radius))
            outer_circle = Circle(Element(center=center, radius=outer_radius))
        else:
            inner_circle = Arc(Element(center=center, radius=inner_radius,
                                       start_angle=startangle*180/np.pi,
                                       end_angle=endangle*180/np.pi))
            outer_circle = Arc(Element(center=center, radius=outer_radius,
                                       start_angle=startangle*180/np.pi,
                                       end_angle=endangle*180/np.pi))

        new_elements = []
        for e in self.elements(Shape):
            if isinstance(e, Line):
                new_elements += self.copy_line(center, radius, startangle, endangle,
                                               start_line, end_line, inner_circle, outer_circle, e)
                    
            elif isinstance(e, Arc):
                new_elements += self.copy_arc(center, radius, startangle, endangle,
                                              start_line, end_line, inner_circle, outer_circle, e)

            elif isinstance(e, Circle):
                new_elements += self.copy_circle(center, radius, startangle, endangle,
                                                 start_line, end_line, inner_circle, outer_circle, e)

        return Geometry(new_elements, adjust=False)   

    def is_new_angle(self, alpha_list, alpha):
        for a in alpha_list:
            if np.isclose(a, alpha):
                return False
        return True

    def search_symmetrie_lines(self, center, radius):
        alpha_list = []
        for e in self.elements(Line):
            alpha_p1 = alpha_line(center, e.p1)
            alpha_p2 = alpha_line(center, e.p2)
            if np.isclose(alpha_p1, alpha_p2):
                if self.is_new_angle(alpha_list, alpha_p1):
                    alpha_list.append(alpha_p1)
        return alpha_list

    def is_symmetrie_line(self, line):
        points = []
        for e in self.elements(Shape):
            points += e.intersect_line(line, self.pickdist, False )
        if len(points) > 0:
            return False
        return True
           
    def find_symmetrie(self, center, radius):
        self.symlines = []
        for alpha in self.search_symmetrie_lines(center, radius):
            p = point(center, radius+5, alpha)
            line = Line(Element(start=center, end=p))
            if self.is_symmetrie_line(line):
                self.add_schnitt(line)
                self.symlines.append(alpha)
        
    def __str__(self):
        return "Nodes {}\nConnected {}".format(
            len(self.g.nodes()),
            nx.number_connected_components(self.g))

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self.__dict__ == other.__dict__)

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)
    
    def __hash__(self):
        """Override the default hash behavior
        (that returns the id or the object)"""
        return hash(tuple(sorted(self.__dict__.items())))
    
    def intersect(self, line, ndt=1.0):
        """returns intersection points of edges with line"""
        path = np.hstack([e[2]['object'].range(ndt)
                          for e in self.g.edges(data=True)])
        logger.debug("path len %d", len(path[0]))
        if np.isclose(line.dx(), 0, self.pickdist):
            idx = np.argwhere(np.isclose(path[0],
                                         line.p1[0], atol=ndt/2)).reshape(-1)
        else:
            idx = np.argwhere(np.isclose(path[1],
                                         line(path[0]),
                                         atol=ndt/2)).reshape(-1)
        return np.asarray(path[0])[idx], np.asarray(path[1])[idx]

    def add_arc(self, center, radius, start_angle, end_angle, ndt=1.0):
        """adds an arc without splitting edges"""
        arc = Arc(Element(center=center, radius=radius,
                  start_angle=start_angle*180/np.pi,
                  end_angle=end_angle*180/np.pi))
        n = self.find_nodes(arc.start(), arc.end())
        self.g.add_edge(n[0], n[1], object=arc)

    def add_line(self, p0, p1, split=True):
        """adds a line from p0 to p1 and splits edges at its
        intersection points"""
        if not split:
            n = self.find_nodes(p0, p1)
            self.g.add_edge(n[0], n[1],
                            object=Line(Element(start=n[0],
                                                end=n[1])))
            return
        line = Line(Element(start=p0, end=p1))
        intersect = []
        edges = []
        # collect all intersecting edges with line
        for e in self.g.edges(data=True):
            #logger.info("range %s", e[2]['object'])
            for ip in e[2]['object'].intersect(line, self.pickdist,
                                               include_end=True):
                if ip:
                    logger.info("append intersection %s", ip)
                    intersect.append(ip)
                    edges.append(e)

        logger.info("edges %d", len(edges))
        # split and add edges and collect intersection points
        for e in edges:
            newshapes = e[2]['object'].split(line, self.pickdist)
            if newshapes:
                logger.info("new shapes %s", len(newshapes))
                self.g.remove_edge(e[0], e[1])
                for s in newshapes:
                    if not np.isclose(s.start(), s.end(),
                                      self.pickdist).all():
                        logger.info("adding edge %s -- %s",
                                    s.start(), s.end())
                        n = self.find_nodes(s.start(), s.end())
                        self.g.add_edge(n[0], n[1], object=s)
                        
        # split all circles and collect intersection points
        for c in self.circles():
            newshapes = c.split(line, self.pickdist)
            if newshapes:
                logger.info("new circle shapes %s", len(newshapes))
                self.g.remove_node(c.center)
                for s in newshapes:
                    if not np.isclose(s.start(), s.end(),
                                      self.pickdist).all():
                        logger.info("adding edge %s -- %s",
                                    s.start(), s.end())
                        n = self.find_nodes(s.start(), s.end())
                        self.g.add_edge(n[0], n[1], object=s)

                    
        # add the line from p0 to p1 
        logger.info("intersections %s", set(intersect))
        for p in sorted(set(intersect), reverse=np.greater(p0, p1).any()):
            if not np.isclose(p0, p, self.pickdist).all():
                logger.info("  %s -- %s", p0, p)
                line = Line(Element(start=p0, end=p))
                n = self.find_nodes(line.start(), line.end())
                self.g.add_edge(n[0], n[1], object=line)
            p0 = p
        if not np.isclose(p0, p1, self.pickdist).all():
            line = Line(Element(start=p0, end=p1))
            n = self.find_nodes(line.start(), line.end())
            self.g.add_edge(n[0], n[1], object=line)

    def set_name(self, name):
        self._name = name
        
    def get_name(self):
        return self._name
        
    def min(self):
        return min(self.corners)

    def max(self):
        return max(self.corners)

    def argmax(self):
        a = np.asarray(self.corners).T
        return np.max(np.arctan2(a[1], a[0]))

    def minmax(self):
        """ Die Funktion bestimmt das Minimum und Maximum auf der x- und der
            y-Achse über alle Shape-Objekte in der Geometrie
            (return [<min-x>, <max-x>, <min-y>, <max-y>])
        """
        mm = [99999,-99999,99999,-99999]
            
        for e in self.elements(Shape):
            n = e.minmax()
            mm[0] = min(mm[0], n[0])
            mm[1] = max(mm[1], n[1])
            mm[2] = min(mm[2], n[2])
            mm[3] = max(mm[3], n[3])
            
        return mm

    def add_schnitt(self, e):
        self.schnitt.append(e)

    def clear_schnitt(self):
        self.schnitt = []
        
    def render_schnitt(self, renderer):
        for e in self.schnitt:
            e.render(renderer, 'darkred')

    def render_airgaps(self, renderer):
        for g in self.airgaps:
            if isinstance(g, Arc):
                renderer.arc(g.startangle, g.endangle, g.center, g.radius, color='red')
            elif isinstance(g, Circle):
                renderer.circle(g.center, g.radius, color='red')
            elif isinstance(g, Line):
                renderer.line(g.p1, g.p2, color='red')
            elif isinstance(g, Point):
                renderer.point(g.p1, 'ro')

    def check_hull(self, center, radius, x, y, pickdist):
        node_count = 0
        miss_count = 0
        for h in convex_hull(self.virtual_nodes()):
            dist = distance(center, h)
            node_count+=1
            if not np.isclose(dist, radius, pickdist):
                if x != None:
                    if np.isclose(x, h[0], pickdist, 1e-04):
                        continue
                if y != None:
                    if np.isclose(y, h[1], pickdist, 1e-04):
                        continue
                miss_count+=1
                
        return miss_count == 0
        
    def get_motor(self):
        mm = self.minmax()
        height = mm[3]-mm[2]
        width = mm[1]-mm[0]
        print("\nGet Motor mit minmax={}".format(mm))
        print("   Motor hoch={}, breit={}".format(height, width))
        
        c = []
        r = 0.0
        if np.isclose(height, width, self.pickdist):
            r = width/2
            c = [mm[1]-r, mm[3]-r]
            logger.info("check for full motor")
            if self.check_hull(c, r, None, None, self.pickdist):
                return Motor(c, r, 0.0, 0.0)
                
            logger.info("check for quarter motor")
            r = width
            c = [mm[0], mm[2]]
            if self.check_hull(c, r, mm[0], mm[2], self.pickdist):
                return Motor(c, r, 0.0, np.pi/2)
            
        elif np.isclose(width, height*2, self.pickdist):
            r = width/2
            c = [mm[1]-height, mm[2]]
            logger.info("check for half motor")
            if self.check_hull(c, r, None, mm[2], self.pickdist):
                return Motor(c, r, 0.0, np.pi)

            c = [mm[1]-height, mm[3]]
            if self.check_hull(c, r, None, mm[3], self.pickdist):
                return Motor(c, r, np.pi, 0.0)

        elif np.isclose(width*2, height, self.pickdist):
            r = width
            c = [mm[0], mm[1]-width]
            logger.info("check for half motor")
            c = [mm[1], mm[3]-width]
            if self.check_hull(c, r, mm[1], None, self.pickdist):
                return Motor(c, r, np.pi/2.0, -np.pi/2.0)

            c = [mm[0], mm[3]-width]
            if self.check_hull(c, r, mm[0], None, self.pickdist):
                return Motor(c, r, -np.pi/2.0, np.pi/2.0)

        motor = self.get_motor_part(mm)
        if motor != None:
            return motor
            
        return Motor([0.0, 0.0], 0.0, 0.0, 0.0 )

    def is_new_center(self, center_list, center, pickdist):
        for c in center_list:
            if points_are_close(c[0], center, pickdist):
                c[1][0] += 1
                return False
        return True
        
    def get_motor_part(self, mm):
        print("get_motor_part")
        center_list = []
        for e in self.elements(Arc):
            center = [round(e.center[0],3), round(e.center[1],3)]
            if self.is_new_center(center_list, center, self.pickdist):
                center_list.append((center, [1]))

        center = []
        count = 0
        unique = False
        for c in center_list:
            if c[1][0] == count:
                unique = False
            elif c[1][0] > count:
                count = c[1][0]
                unique = True
                center = c[0]
                
        if not unique:
            # Wir finden keine Arc-Objekte, welche uns einen Hinweis auf den
            # Center geben können. Wir versuchen in der Verzweiflung mit
            # x(min) und y(min)
            print("what's that?")
            center = [round(mm[0], 4), round(mm[2], 4)]

        print("Center is {}".format(center))                
        
        min_radius = 99999
        max_radius = 0
        startangle = 999.0
        endangle = -999.0
        
        for h in convex_hull(self.virtual_nodes()):
            angle = alpha_line(center, [round(h[0], 4), round(h[1], 4)])
            if angle < 0.0:
                print("Seltsamer Punkt {}".format(h))
            startangle = min(startangle, angle)
            endangle = max(endangle, angle)
            
            dist = distance(center, h)
            min_radius = min(min_radius, dist)
            max_radius = max(max_radius, dist)

        center_left = round(center[0] - mm[0], 4)
        center_right = round(mm[1] - center[0], 4)
        center_down = round(center[1] - mm[2], 4)
        center_up = round(mm[3] - center[1], 4)
        
        min_r = min(center_left, center_right, center_up, center_down)
        max_r = max(center_left, center_right, center_up, center_down)
        
        print("   left={}, right={}, up={}, down={}"
            .format(center_left, center_right, center_up, center_down))        
        print("   startangle={}, endangle={}".format(startangle, endangle))
        print("   min_radius={}, max_radius={}".format(min_radius, max_radius))


        if less_equal(center_left, 0.0) or less_equal(center_right, 0.0) or \
           less_equal(center_up, 0.0) or less_equal(center_down, 0.0):
            if not np.isclose(center_left, center_right):
                x = part_of_circle(startangle, endangle)
                print("Teiler {}, start={}, end={}".format(x, startangle, endangle))

                if x > 2:
                    print("Teiler {} => ok".format(x))
                    return Motor([round(center[0], 8), round(center[1], 8)], max_radius, startangle, endangle)
        
        if min_radius >= max_radius*0.9 and min_r >= max_r*0.9:
            # Mit 10 % Abweichungen gehen wir noch von einem ganzen Motor aus.
            return Motor([round(center[0], 8), round(center[1],8)], max(max_radius, max_r), 0.0, 0.0 )

        if np.isclose(center_down, 0.0):
            min_r = min(center_left, center_right, center_up)
            if min_r >= max_r*0.9:
                # Vermutlich ein halber Motor
                return Motor([round(center[0], 8), round(center[1], 8)], max(max_radius, max_r), 0.0, np.pi )
                
            if np.isclose(center_left, 0.0):
                min_r = min(center_right, center_up)
                if min_r >= max_r*0.9:
                    # Vermutlich ein viertel Motor
                    return Motor([round(center[0], 8), round(center[1],8)], max(max_radius, max_r), 0.0, np.pi/2 )

        # Die Varianten von halben und viertel Motoren mit anderer Drehung sparen
        # wir uns für später auf.

        if np.isclose(center_left, center_right):
            # Der x-Wert des Center scheint zu stimmen. Eine Symmetrie an der y-Achse
            # ist möglich. Dem y-Wert ist nicht zu trauen!
            if center_down > 0.0:
                # Der Center ist zu hoch! Wir müssen neu rechnen.
                nodes = [n for n in convex_hull(self.g.nodes()) if n[0] > center[0]]
                assert(nodes)
                p = nodes[0]
                for n in nodes[1:]:
                    if n[1] < p[1]:
                        p = n
                    
                m_min = 99999.0
                for n in nodes:
                    m = line_m(p, n)
                    if m != None and m > 0.0:
                        #print("m = {}, p={}, n={}".format(m, p, n))
                        m_min = min(m_min, m)
                    
                y = line_n([p[0]-center[0], p[1]], m_min)
                center[1] = y
                angle = alpha_line(center,p)

        return Motor([round(center[0], 8), round(center[1], 8)], max_radius, angle, np.pi - angle)
    
    def is_new_radius(self, radius_list, radius):
        for r,d in radius_list:
            if np.isclose(r, radius):
                return False
        return True

    def is_airgap(self, center, radius, startangle, endangle, circle):
        """ Die Funktion untersucht, ob sich der Parameter circle in einem
            Luftspalt befindet.
        """
        ok=True
        for e in self.elements(Shape):
            for p in e.intersect_circle(circle, 0.0, True):
                alpha_p = alpha_line(center, p)
                if not (np.isclose(alpha_p, startangle) or
                        np.isclose(alpha_p, endangle)):
                    self.airgaps.append(Point(p))
                    ok=False                          
        return ok

    def is_border_line(self, center, startangle, endangle, e):
        if isinstance(e, Line):
            angle_p1 = alpha_line(center, e.p1)
            if np.isclose(startangle, angle_p1):
                angle_p2 = alpha_line(center, e.p2)
                return np.isclose(startangle, angle_p2)
            elif np.isclose(endangle, angle_p1):
                angle_p2 = alpha_line(center, e.p2)
                return np.isclose(endangle, angle_p2)
        return False
        
    def detect_airgaps(self, center, startangle, endangle):
        """ Die Funktion sucht Luftspalt-Kandidaten und liefert eine Liste
            von Möglichkeiten mit jeweils einem minimalen und einem maximalen
            Radius als Begrenzung des Luftspalts.
        """
        print("detect_airgaps: center={}".format(center))
        
        gaplist = []
        for e in self.elements(Shape):
            if not self.is_border_line(center, startangle, endangle, e):
                gaplist += [e.minmax_from_center(center)]
        gaplist.sort()

        airgaps = []
        dist_max = 0.0
        
        min_radius = gaplist[0][0]
        max_radius = gaplist[len(gaplist)-1][1]
        dist_max = min_radius
        
        for g in gaplist:
            if not less_equal(g[0], dist_max):
                if not (np.isclose(min_radius, dist_max, 1e-2, 1.0) or \
                        np.isclose(g[0], max_radius, 1e-2, 1.0)):
                    airgaps.append((dist_max, g[0]))
            dist_max = max(dist_max ,g[1])
            
        return airgaps