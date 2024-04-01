import logging
import numpy as np
from .shape import Circle, Arc, Line, Element
from .functions import distance

logger = logging.getLogger(__name__)

def polylines(entity, lf, rf, xoff=0.0, yoff=0.0, rotation=0.0):
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
        except Exception:
            if not entity.is_closed:
                break
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
            yield Arc(Element(center=(pc[0], pc[1]),
                              radius=np.abs(r),
                              start_angle=sa/rf,
                              end_angle=ea/rf),
                      lf, rf,
                      xoff=xoff, yoff=yoff,
                      rotation=rotation)
        else:
            yield Line(Element(start=p1, end=p2), lf,
                       xoff=xoff, yoff=yoff,
                       rotation=rotation)
        i += 1


def lw_polyline(entity, lf, xoff=0.0, yoff=0.0, rotation=0.0):
    """returns a collection of bulged vertices
    http://www.afralisp.net/archive/lisp/Bulges1.htm
    """
    if isinstance(entity.points, list):
        points = [(lf*p[0], lf*p[1]) for p in entity.points]
    else:
        points = [(lf*p[0], lf*p[1]) for p in entity.points()]

    if points:
        p1 = points[0]
        for p2 in points[1:]:
            yield Line(Element(start=p1, end=p2), lf,
                       xoff=xoff, yoff=yoff,
                       rotation=rotation)
            p1 = p2
    if entity.is_closed:
        yield Line(Element(start=p1, end=points[0]), lf,
                   xoff=xoff, yoff=yoff,
                   rotation=rotation)


def ellipse(entity, lf, xoff=0.0, yoff=0.0, rotation=0.0):
    w = np.linalg.norm(entity.major_axis) * 2
    h = entity.ratio * w
    theta = np.arctan2(entity.major_axis[1], entity.major_axis[0])
    start_angle = entity.start_param
    end_angle = entity.end_param
    if end_angle < start_angle:
        end_angle += 2*np.pi
    alfa = np.linspace(start_angle, end_angle, 20)
    x = 0.5 * w * np.cos(alfa)
    y = 0.5 * h * np.sin(alfa)
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])
    x, y = np.dot(R, [x, y])
    x += entity.center[0]
    y += entity.center[1]
    points = np.array((x, y)).T
    p1 = points[0]
    for p2 in points[1:]:
        yield Line(Element(start=p1, end=p2), lf,
                   xoff=xoff, yoff=yoff,
                   rotation=rotation)
        p1 = p2


def spline(entity, lf, min_dist=0.001, xoff=0.0, yoff=0.0, rotation=0.0):
    if False:
        yield Line(Element(start=entity.control_points[0],
                           end=entity.control_points[-1]), lf,
                   xoff=xoff, yoff=yoff,
                   rotation=rotation)
        return

    if False:
        p_prev = None
        for p in entity.control_points:
            if p_prev:
                yield Line(Element(start=p_prev, end=p), lf,
                           xoff=xoff, yoff=yoff,
                           rotation=rotation)
            p_prev = p
        return

    points_between = entity.control_points[1:-1]
    p1 = entity.control_points[0]
    pe = entity.control_points[-1]
    for p2 in points_between:
        dist_12 = distance(p1, p2)
        dist_2e = distance(p2, pe)
        if dist_2e < min_dist:
            logger.debug("SPLINE: ignor small end-distance %s", dist_2e)
            yield Line(Element(start=p1, end=pe), lf,
                       xoff=xoff, yoff=yoff,
                       rotation=rotation)
            return

        if dist_12 > min_dist:
            yield Line(Element(start=p1, end=p2), lf,
                       xoff=xoff, yoff=yoff,
                       rotation=rotation)
            p1 = p2
        else:
            logger.debug("SPLINE: ignor small distance %s", dist_12)

    yield Line(Element(start=p1, end=pe), lf,
               xoff=xoff, yoff=yoff,
               rotation=rotation)


def face3d(entity, lf):
    logger.info("FACE3D: Points=%s", entity.points)
    for i in range(len(entity.points)-1):
        if not entity.is_edge_invisible(i):
            ip = i+1 if i < 4 else 0
            yield Line(Element(start=(entity.points[i][1],
                                      entity.points[i][2]),
                               end=(entity.points[ip][1],
                                    entity.points[ip][2])))


def insert_block(dwg, insert_entity, lf, rf, block, min_dist=0.001):
    logger.debug('Insert %s entities from block %s',
                 len(block),
                 insert_entity.name)
    logger.debug('Insert = %s', insert_entity.insert)
    logger.debug('Rotation = %s', insert_entity.rotation)
    logger.debug('Scale = %s', insert_entity.scale)
    logger.debug('Rows = %s', insert_entity.row_count)
    logger.debug('Cols = %s', insert_entity.col_count)
    logger.debug('Row spacing = %s', insert_entity.row_spacing)
    logger.debug('Col spacing = %s', insert_entity.col_spacing)

    if insert_entity.insert != (0.0, 0.0, 0.0):
        logger.debug('Different Location in Insert')

    xoff = insert_entity.insert[0]
    yoff = insert_entity.insert[1]

    scale = (round(insert_entity.scale[0], 8),
             round(insert_entity.scale[1], 8),
             round(insert_entity.scale[2], 8))
    if not (scale == (1.0, 1.0, 1.0) or
            scale == (1.0, 1.0, 0.0)):
        logger.error('Block scaling in Insert not supported')
        logger.error('  scale = {}'.format(scale))
        return

    if(insert_entity.row_count > 1 or
       insert_entity.col_count > 1 or
       insert_entity.row_spacing > 0 or
       insert_entity.col_spacing > 0):
        logger.error('Multi Block references in Insert not supported')
        return

    for e in block:
        if e.dxftype == 'ARC':
            logger.debug("Block Arc")
            yield Arc(e, lf, rf,
                      xoff=xoff, yoff=yoff,
                      rotation=insert_entity.rotation)
        elif e.dxftype == 'CIRCLE':
            logger.debug("Block Circle %s, Radius %f", e.center[:2], e.radius)
            yield Circle(e, lf,
                         xoff=xoff, yoff=yoff,
                         rotation=insert_entity.rotation)
        elif e.dxftype == 'LINE':
            logger.debug("Block Line")
            yield Line(e, lf,
                       xoff=xoff, yoff=yoff,
                       rotation=insert_entity.rotation)
        elif e.dxftype == 'POLYLINE':
            logger.debug("Block Polyline")
            for p in polylines(e, lf, rf,
                               xoff=xoff, yoff=yoff,
                               rotation=insert_entity.rotation):
                yield p
        elif e.dxftype == 'LWPOLYLINE':
            logger.debug("Block lwPolyline")
            for p in lw_polyline(e, lf,
                                 xoff=xoff, yoff=yoff,
                                 rotation=insert_entity.rotation):
                yield p
        elif e.dxftype == 'SPLINE':
            logger.debug("Block spline")
            for l in spline(e, lf,
                            min_dist=min_dist,
                            xoff=xoff, yoff=yoff,
                            rotation=insert_entity.rotation):
                yield l
        elif e.dxftype == 'INSERT':
            logger.debug("Nested Insert of Block %s", e.name)
            block = dwg.blocks[e.name]
            for l in insert_block(dwg, e, lf, rf, block, min_dist=min_dist):
                yield l

        else:
            logger.warn("unknown type %s in block %s",
                        e.dxftype, insert_entity.name)


def dxfshapes0(dxffile, mindist=0.01, layers=[]):
    """returns a collection of dxf entities (ezdxf)"""
    import ezdxf
    dwg = ezdxf.readfile(dxffile)
    id = 0
    # $ACADVER: AC1006 = R10, AC1009 = R11 and R12, AC1012 = R13,
    #   AC1014 = R14 AC1015 = Release 2000/0i/2
    # check units:
    # dwg.header['$ANGDIR'] 1 = Clockwise angles, 0 = Counterclockwise
    # dwg.header['$AUNITS'] 0 Decimal Degrees, 1 Deg/Min/Sec, 2 Grads, 3 Radians
    # dwg.header['$INSUNIT'] 1 = Inches; 2 = Feet; 3 = Miles;
    #   4 = Millimeters; 5 = Centimeters; 6 = Meters
    # dwg.header['$LUNITS']
    for e in dwg.modelspace():
        if e.dxftype() == 'ARC':
            yield Arc(e.dxf)
        elif e.dxftype() == 'CIRCLE':
            logger.debug("Circle %s, Radius %f", e.center[:2], e.radius)
            yield Circle(e.dxf)
        elif e.dxftype() == 'LINE':
            yield Line(e.dxf)
        elif e.dxftype() == 'POLYLINE':
            for p in polylines(e):
                yield p
        elif e.dxftype() == 'SPLINE':
            for l in spline(e, 1.0, in_dist=mindist):
                yield l
        elif e.dxftype() == 'POINT':
            logger.debug("Id %d4: type %s ignored", id, e.dxftype)
        else:
            logger.warning("Id %d4: unknown type %s", id, e.dxftype)
        id += 1


def dxfshapes(dxffile, mindist=0.01, layers=[]):
    """returns a collection of dxf entities (dxfgrabber)"""
    import dxfgrabber
    dwg = dxfgrabber.readfile(dxffile)
    # print("Layers = {}".format(dwg.layers.names()))
    id = 0
    # $ACADVER: AC1006 = R10, AC1009 = R11 and R12, AC1012 = R13,
    #   AC1014 = R14 AC1015 = Release 2000/0i/2
    # check units:
    # dwg.header['$ANGDIR'] 1 = Clockwise angles, 0 = Counterclockwise
    # dwg.header['$AUNITS'] Decimal Degrees, Deg/Min/Sec, Grads, Radians
    # dwg.header['$INSUNIT'] 1 = Inches; 2 = Feet; 3 = Miles;
    #   4 = Millimeters; 5 = Centimeters; 6 = Meters
    # dwg.header['$LUNITS']
    lf = 1
    if dwg.header.get('$LUNITS', 0) == 1:
        # conv = [1, 2.54e-2, 10.12, 633.0, 1e-3, 1e-2, 1]
        lf = 2.54e3

    rf = np.pi/180
    if dwg.header.get('$AUNITS', 0) == 4:
        rf = 1

    for e in dwg.modelspace():
        if not layers or e.layer in layers:
            if e.dxftype == 'ARC':
                yield Arc(e, lf, rf)
            elif e.dxftype == 'CIRCLE':
                logger.debug("Circle %s, Radius %f", e.center[:2], e.radius)
                yield Circle(e, lf)
            elif e.dxftype == 'LINE':
                yield Line(e, lf)
            elif e.dxftype == 'POLYLINE':
                for p in polylines(e, lf, rf):
                    yield p
            elif e.dxftype == 'LWPOLYLINE':
                for p in lw_polyline(e, lf):
                    yield p
            elif e.dxftype == 'SPLINE':
                for l in spline(e, lf, min_dist=mindist):
                    yield l
            elif e.dxftype == 'INSERT':
                logger.debug("Insert of Block %s", e.name)
                block = dwg.blocks[e.name]
                for l in insert_block(dwg, e, lf, rf, block, min_dist=mindist):
                    yield l
            elif e.dxftype == 'ELLIPSE':
                for l in ellipse(e, lf):
                    yield l
                #w = np.linalg.norm(e.major_axis) * 2
                #h = e.ratio * w
                #rtheta = np.arctan2(e.major_axis[1], e.major_axis[0])
                #angle = rtheta*180/np.pi
                #start_angle = e.start_param*180/np.pi + angle
                #end_angle = e.end_param*180/np.pi + angle
                #if end_angle < start_angle:
                #    end_angle += 360
                #arc = Arc(Element(center=e.center,
                #                  radius=w/2,
                #                  start_angle=start_angle,
                #                  end_angle=end_angle,
                #                  width=w,
                #                  height=h,
                #                  rtheta=rtheta,
                #                  start_param=e.start_param,
                #                  end_param=e.end_param))
                #yield arc

            elif e.dxftype == 'POINT':
                logger.debug("Id %d4: type %s ignored", id, e.dxftype)
            elif e.dxftype == '3DFACE':
                logger.warning(
                    "Id %d4: type %s not implemented", id, e.dxftype)
                # for l in face3d(e, lf):
                #     yield l
            else:
                logger.warning("Id %d4: unknown type %s", id, e.dxftype)
            id += 1
