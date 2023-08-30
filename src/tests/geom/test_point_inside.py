import femagtools.dxfsl.geom as g

import logging


area = [[
    g.Arc(g.Element(center=(93.000391, 104.825783),
                    radius=90.000363,
                    start_angle=-113.674935,
                    end_angle=-102.876410)),
    g.Arc(g.Element(center=(72.676972, 15.919097),
                    radius=1.199853,
                    start_angle=12.363429,
                    end_angle=77.142632)),
    g.Arc(g.Element(center=(0.000000, 0.000000),
                    radius=75.599833,
                    start_angle=10.725943,
                    end_angle=12.355028)),
    g.Arc(g.Element(center=(73.542349, 13.930019),
                    radius=0.749743,
                    start_angle=-99.701978,
                    end_angle=10.759233)),
    g.Arc(g.Element(center=(90.226618, 111.767627),
                    radius=99.999841,
                    start_angle=-110.995630,
                    end_angle=-99.677752)),
    g.Arc(g.Element(center=(75.894867, 74.423271),
                    radius=59.999926,
                    start_angle=-134.554066,
                    end_angle=-110.995752)),
    g.Arc(g.Element(center=(34.501683, 32.380386),
                    radius=0.999352,
                    start_angle=135.002787,
                    end_angle=-134.566333)),
    g.Line(g.Element(start=(36.133000, 35.426000),
                     end=(33.795000, 33.087000))),
    g.Arc(g.Element(center=(36.840535, 34.719239),
                    radius=1.000084,
                    start_angle=45.700702,
                    end_angle=135.031345)),
    g.Arc(g.Element(center=(82.961666, 81.929844),
                    radius=64.999916,
                    start_angle=-134.331700,
                    end_angle=-113.674891))
],
    # deeg example
    [g.Arc(g.Element(center=(2., 6.),
                     radius=1.,
                     start_angle=0.,
                     end_angle=180.)),
     g.Arc(g.Element(center=(6., 2.),
                     radius=1.,
                     start_angle=-90.,
                     end_angle=90.)),
     g.Arc(g.Element(center=(4., 4.),
                     radius=1.,
                     start_angle=180.,
                     end_angle=270.)),
     g.Arc(g.Element(center=(4., 4.),
                     radius=3.,
                     start_angle=180.,
                     end_angle=-90.)),
     g.Line(g.Element(start=(1., 6.), end=(1., 4.0))),
     g.Line(g.Element(start=(3., 4.), end=(3., 6.0))),
     g.Line(g.Element(start=(6., 1.), end=(4., 1.0))),
     g.Line(g.Element(start=(4., 3.), end=(6., 3.0)))]
]


def check_point_inside(n):
    geom = g.Geometry(area[n])
    geom.center = (0.0, 0.0)

    pts = []
    for a in geom.list_of_areas():
        pt = a.get_point_inside(geom)
        assert(pt)
        pts.append(pt)

    assert(pts)
#    p = dr.PlotRenderer()
#    p.render_elements(geom, g.Shape, neighbors=True, points=pts)


def test_0_point_inside():
    check_point_inside(0)


def test_1_point_inside():
    check_point_inside(1)

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    test_1_point_inside()
