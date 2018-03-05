import pytest
import numpy as np
from femagtools import isa7

def test_read():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)
    
    assert isa.NUM_PNT == 39
    assert isa.NUM_LIN == 39
    assert isa.NUM_SPEL_NDCH == 112
    assert isa.POINT_ISA_PT_VALID[0] == True
    assert type(isa.POINT_ISA_PT_VALID[0]) == bool
    assert isa.POINT_ISA_POINT_REC_PT_CO_X[0] == pytest.approx(0.008339,
                                                               abs=1e-5)
    
def test_objects():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)

    assert len(isa.points) == 39
    assert len(isa.lines) == 39
    assert len(isa.nodes) == 1729
    assert len(isa.nodechains) == 98
    assert len(isa.elements) == 822
    assert len(isa.superelements) == 3

def test_points():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)
    p = isa.points[0].x, isa.points[0].y
    p_xy = isa.points[0].xy

    assert p == p_xy == pytest.approx((0.008339, 0.011702), abs=1e-5)

def test_lines():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)
    ln = isa.lines[0]

    assert type(ln.p1) == isa7.Point
    assert type(ln.p2) == isa7.Point
    assert ln.p1.xy == pytest.approx((0.008339, 0.011702), abs=1e-5)
    assert ln.p2.xy == pytest.approx((0.007388, 0.007190), abs=1e-5)

def test_nodes():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)
    nd = isa.nodes[-1]

    assert type(nd) == isa7.Node
    assert nd.xy == (nd.x, nd.y) == pytest.approx(
        (0.008787, 0.007736), abs=1e-5)

def test_nodechains():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)
    nc = isa.nodechains[0]

    assert type(nc) == isa7.NodeChain
    assert type(nc.nodes[0]) == isa7.Node

def test_nodechains_reverse():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)

    for nc in isa.nodechains:
        assert nc.reverse().nodes == nc.nodes[::-1]

def test_elements():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)
    el = isa.elements[0]

    assert type(el) == isa7.Element
    for v in el.vertices:
        assert type(v) == isa7.Node
    assert el.vertices[0].xy == pytest.approx((0.002202, 0.003197),
                                              abs=1e-5)

def test_superelements():
    filename = 'tests/data/minimal.ISA7'
    isa = isa7.read(filename)
    se = isa.superelements[0]

    assert type(se) == isa7.SuperElement
    for el in se.elements:
        assert type(el) == isa7.Element
    for nc in se.nodechains:
        assert type(nc) == isa7.NodeChain
        assert nc.key == se.nc_keys[se.nodechains.index(nc)]

def test_subregions():
    filename = 'tests/data/test_disp_stat.ISA7'
    isa = isa7.read(filename)
    sr = isa.subregions[0]

    assert type(sr) == isa7.SubRegion
    assert sr.name == "Stat"
    for se in sr.superelements:
        assert type(se) == isa7.SuperElement
    for nc in sr.nodechains:
        assert type(nc) == isa7.NodeChain

def test_windings():
    filename = 'tests/data/test_disp_stat.ISA7'
    isa = isa7.read(filename)
    wd = isa.windings[0]
    
    assert type(wd) == isa7.Winding
    assert wd.name == "Stra"
    for sr in wd.subregions:
        assert type(sr) == isa7.SubRegion
    assert wd.num_turns == 100
