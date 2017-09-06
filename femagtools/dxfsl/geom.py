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

class Corner(object):
    def __init__(self, p):
        self.__p = p
        self.__keep = False

    def point(self):
        return self.__p
        
    def is_equal(self, p):
        return np.isclose(self.__p[0], p[0]) and np.isclose(self.__p[1], p[1])
        
    def length(self):
        return np.sqrt(self.__p[0]**2 + self.__p[1]**2)

    def same_angle(self, p):
        return self.angle() == Corner(p).angle()

    def angle(self):
        arg = np.arctan2(self.__p[1], self.__p[0]) * 180. / np.pi
        while arg < 0.0:
            arg = 360. + arg
        return round(arg, 6)

    def set_keep_node(self):
        self.__keep = True

    def keep_node(self):
        return self.__keep

    def __lt__(self, c):
        if self.angle() < c.angle():
            return 1
        if self.angle() > c.angle():
            return 0
                
        if self.length() < c.length():
            return 1
        else:
            return 0
        
    def __str__(self):
        return "Corner: ({},{}) -len={} -angle={}".format(self.__p[0],
                                                          self.__p[1],
                                                          self.length())


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

    def get_p1(self):
        return Point(self.p1[0], self.p1[1])
        
    def get_p2(self):
        return Point(self.p2[0], self.p2[1])
        
    def get_other_point(self, pt):
        if pt.x == self.p1[0] and pt.y == self.p1[1]:
            return Point(self.p2[0], self.p2[1])

        if pt.x == self.p2[0] and pt.y == self.p2[1]:
            return Point(self.p1[0], self.p1[1])
        
        return None
        
    def get_center(self):
        p = self.center_of_connection()
        return Point(p[0], p[1])

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
    
    def __str__(self):
        return " {}/{}".format(self.p1, self.p2)
    

class Circle(Shape):
    """a circle with center and radius"""
    def __init__(self, e):
        self.center = e.center[:2]
        self.radius = e.radius
        self.p1 = self.center[0]-self.radius, self.center[1]
        self.p2 = self.center[0]+self.radius, self.center[1]
        
    def render(self, renderer):
        renderer.circle(self.center, self.radius)

    def move(self, dist):
        super(Circle, self).move(dist)
        self.center = self.center[0]+dist[0], self.center[1]+dist[1]
        
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

    def split(self, x, pickdist):
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
                            
    def render(self, renderer):
            renderer.arc(self.startangle, self.endangle,
                         self.center, self.radius)
        
        #num = int(self.length()/ndt) if ndt > 0 else 0
        #return u"nc_circle_m({}, {}, {}, {}, {}, {}, {})\n".format(
        #    self.p1[0], self.p1[1], self.p2[0], self.p2[1],
        #    self.center[0], self.center[1], num)

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
    
    def split(self, x, pickdist):
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
    
    def is_inside(self, alpha, pickdist, include_end=False):
        """returns True if alpha is between start and end angle"""
        logger.debug("Arc alpha %f in (%f, %f)", alpha,
                     self.startangle, self.endangle)
        a0 = self.startangle
        da = pickdist/self.radius
        a1 = self.endangle
        if (a0+da < alpha < a1-da):
            return True
        if (a0+da < alpha+2*np.pi < a1-da):
            return True
        if include_end:
            if (np.less_equal(a0-da, alpha) and
                np.less_equal(alpha, a1+da)):
                return True
            alpha += 2*np.pi
            return (np.less_equal(a0-da, alpha) and
                    np.less_equal(alpha, a1+da))

        return False

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

    def __str__(self):
        return "Arc c {} r {} {} -- {} {}".format(self.center,
                                                  self.radius,
                                                  self.startangle*180/np.pi,
                                                  self.endangle*180/np.pi,
                                                  self.center_of_connection())
    
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


class Line(Shape):
    """straight connection between start and end point"""
    def __init__(self, e):
        self.p1 = e.start[0], e.start[1]
        self.p2 = e.end[0], e.end[1]

    def render(self, renderer):
        renderer.line(self.p1, self.p2)

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
        
    def split(self, line, pickdist):
        """split this line at intersection point with line"""
        point = self.intersect(line, pickdist)
        if point and self.is_inside(point[0], pickdist):
            logger.info("Line intersection %s", point)
            return (Line(Element(start=self.p1, end=point[0])),
                    Line(Element(start=point[0], end=self.p2)))
        return []
    
    def is_inside(self, point, pickdist, include_end=False):
        """returns True if point is between start and end point"""
        logger.debug("%s in (%s, %s)", point, self.p1, self.p2)
        if np.isclose(self.dx(), 0):
            if point[0] - self.p1[0] <= pickdist and \
               point[1] - self.ymin() > pickdist and \
               self.ymax() - point[1] > pickdist:
                return True
              
        elif np.isclose(self.dy(), 0):
            if point[1] - self.p1[1] <= pickdist and \
               point[0] - self.xmin() > pickdist and \
               self.xmax() - point[0] > pickdist:
                return True
                
        elif la.norm(np.asarray(point) - self.p1) > pickdist and \
             la.norm(np.asarray(self.p2) - point) > pickdist and \
             (np.asarray(point) > (self.xmin(), self.ymin())).all() and \
             (np.asarray(point) < (self.xmax(), self.ymax())).all():
            return True
        
        if include_end:
            if la.norm(np.asarray(self.p1) - point) <= pickdist or \
               la.norm(np.asarray(self.p2) - point) <= pickdist:
                return True
        return False
    
    def __str__(self):
        return "Line {} -- {} center {}".format(self.p1,
                                                self.p2,
                                                self.center_of_connection())

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


class Geometry(object):
    """collection of connected shapes"""
    def __init__(self, elements=[], pickdist=1e-3, adjust=True):
        self._name = ''
        self.mirror_axis = None
        self.corners = []
        self.g = nx.Graph()
        self.pickdist = pickdist
        i = 0
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
            
            airgap = self.airgap()
            logger.info("Airgap %s Diameters %s", airgap, self.diameters)
            if airgap:
                self.diameters = (self.diameters[0], round(airgap[0], 3),
                                  round(airgap[1], 3), self.diameters[-1])
            
    def airgap(self):
        """returns airgap diameters if any"""
        nodes = []
        for a in self.area_nodes(incl_bnd=True):
            c = sorted(find_corners(a))
            if len(c) > 2:
                logger.info("Corners %s (%s, %s)", c,
                            self.corners[0], self.corners[-1])
                outer = np.isclose(c, self.corners[-1]).all(axis=1)
                inner = np.isclose(c, self.corners[0]).all(axis=1)
                # check outer:
                if inner.any() and not outer.any():
                    nodes.append(2*round(la.norm(c[-1]), 3))
                elif outer.any() and not inner.any():
                    nodes.append(2*round(la.norm(c[0]), 3))
        try:
            if np.isclose(nodes[0], nodes[1]):
                return []
        except:
            return []
        return sorted(nodes)

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
        rotnodes = np.dot(self.g.nodes(), T)
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

    def fix_hull(self, g):
        """removes corners of g"""
        corners = [Corner(c)
                   for c in self.find_corners(g.nodes())]

        for p1, p2 in g.edges():
            for c in corners:
                if c.is_equal(p1):
                    if not c.same_angle(p2):
                        c.set_keep_node()
                else:
                    if c.is_equal(p2):
                        if c.same_angle(p1):
                            g.remove_edge(p1, p2)
                        else:
                            c.set_keep_node()

        for c in corners:
            if not c.keep_node():
                g.remove_node(c.point())

        corners = [Corner(c) for c in self.find_corners(g.nodes())]
        corners.sort()
            
        # Ein Corner-Objekt mit unmöglichem Winkel
        last_c = Corner([-9, -8])
        for c in corners:
            if last_c.same_angle(c.point()):
                self.add_line(last_c.point(), c.point())
            last_c = c
            
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

    def copy(self, startangle, endangle,
             radius=[0, float("inf")], incl_bnd=True):
        """returns the sector including all nodes and edges
        between start and end angle and radius"""
        clone = Geometry()
        alpha = -startangle
        T = np.array(((np.cos(alpha), -np.sin(alpha)),
                      (np.sin(alpha), np.cos(alpha))))
        for c in self.circles():
            phi = np.arctan2(c.center[1], c.center[0])
            r = la.norm(c.center)
            dphi = np.arctan2(c.radius, r)
            if startangle - dphi <= phi <= endangle + dphi and \
               radius[0] < r + c.radius < radius[1]:
                cc = copy.deepcopy(c).transform(T)
                clone.g.add_node(cc.center, object=cc)
            
        for e in self.g.edges(data=True):
            phi = (np.arctan2(round(e[0][1]/self.pickdist),
                              round(e[0][0]/self.pickdist)),
                   np.arctan2(round(e[1][1]/self.pickdist),
                              round(e[1][0]/self.pickdist)))
            r = la.norm(e[0]), la.norm(e[1])
            dpick = (4*np.arctan2(self.pickdist, r[0]),
                     4*np.arctan2(self.pickdist, r[1]))
            
            logger.debug("r %g %g phi %g, %g (%g, %g)",
                         r[0], r[1], phi[0], phi[1], radius[0], radius[1])
            # check if edge is inside
            if startangle-dpick[0] <= phi[0] <= endangle+dpick[0] and \
               startangle-dpick[1] <= phi[1] <= endangle+dpick[1] and \
               radius[0]-self.pickdist < r[0] < radius[1]+self.pickdist and \
               radius[0]-self.pickdist < r[1] < radius[1]+self.pickdist:
                s = copy.deepcopy(e[2]['object']).transform(T)
                n = clone.find_nodes(s.start(),
                                     s.end())
                logger.debug("%d %s", n)
                add_or_join(clone.g, n[0], n[1], s,
                            clone.pickdist)
                continue

            # check if outside
            if (phi[0] < startangle-dpick[0] and
                phi[1] < startangle-dpick[1]) or \
                (phi[0] > endangle+dpick[0] and
                 phi[1] > endangle+dpick[1]) or \
                 (r[0] < radius[0]-self.pickdist and
                  r[1] < radius[0]-self.pickdist) or \
                 (r[0] > radius[1]+-self.pickdist and
                  r[1] > radius[1]+self.pickdist):
                continue
                
            if (np.asanyarray(r) > radius[1]+self.pickdist).any() or \
               (np.asanyarray(r) < radius[0]-self.pickdist).any():
                continue
            
            # split this edge as one node is outside
            if (np.asarray(phi) > endangle).any():
                r = 1.1*np.sqrt(e[2]['object'].xmax()**2 +
                                e[2]['object'].ymax()**2)
            elif (np.asarray(phi) < startangle).any():
                r = 1.1*np.sqrt(e[2]['object'].xmin()**2 +
                                e[2]['object'].ymin()**2)
            
            if (np.asarray(phi) < startangle).any():
                angle = 0
            else:
                angle = endangle - startangle

            o = copy.deepcopy(e[2]['object']).transform(T)
            line = Line(Element(start=(0, 0),
                                end=(r*np.cos(angle),
                                     r*np.sin(angle))))
            
            for s in o.split(line, self.pickdist):
                cp = s.center_of_connection()
                cphi = np.arctan2(cp[1], cp[0])
                if startangle-dpick[0] < cphi < endangle+dpick[0]:
                    n = clone.find_nodes(s.start(), s.end())
                    logger.info("adding edge %s -- %s", n[0], n[1])
                    clone.g.add_edge(n[0], n[1], object=s)

        clone.corners = sorted(find_corners(clone.g.nodes()))
        clone.set_diameters()
            
        logger.info("Cloned diameters %s", clone.diameters)
        if not clone.diameters:
            return clone
        
        if not incl_bnd:
            #self.remove_corners(clone.g)
            self.fix_hull(clone.g)
            # check boundary if there is more
            logger.info("Check radius %g == %g",
                        self.diameters[-1]/2,
                        la.norm(max(clone.g.nodes())))
            if np.isclose(self.diameters[-1]/2,
                          la.norm(max(clone.g.nodes())),
                          self.pickdist):
                remove_hull(clone.g, self.pickdist)
                
                c = self.find_corners(clone.g.nodes())
                rmin = la.norm(min(c))
                logger.info("Found corners %s rmin %f", c, rmin)
                innerarc = [i for i in c
                            if np.isclose(la.norm(i),
                                          rmin,
                                          5e-3*rmin)]
                rmin = max([la.norm(i) for i in innerarc])
                clone.diameters = (2*rmin, self.diameters[1])
                
                cnodes = get_nodes_of_paths(clone.g, innerarc)
                clone.g.remove_nodes_from(cnodes)
                
#        if len(self.diameters) > 2:
#            clone.diameters = (clone.diameters[0], self.diameters[1],
#                               self.diameters[2], clone.diameters[-1])
        clone.alpha = endangle - startangle
        logger.info("corners %s Diameters %s", clone.corners, clone.diameters)
        dy2 = clone.diameters[0]
        da2 = clone.diameters[-1]
        rot = np.array([[np.cos(clone.alpha), -np.sin(clone.alpha)],
                        [np.sin(clone.alpha), np.cos(clone.alpha)]])
        p0, p1 = np.dot(rot, (da2/2, 0)), np.dot(rot, (dy2/2, 0))
        clone.mirror_axis = p0[0], p0[1], p1[0], p1[1]
        return clone
    
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
