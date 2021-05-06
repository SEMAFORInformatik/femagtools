import pytest
from femagtools import isa7


@pytest.fixture
def isa():
    filename = 'tests/data/minimal.ISA7'
    return isa7.read(filename)


def test_read(isa):
    assert len(isa.points) == 39
    assert len(isa.lines) == 39
    assert isa.points[0].x == pytest.approx(0.008339,
                                            abs=1e-5)
    

def test_objects(isa):
    assert len(isa.points) == 39
    assert len(isa.lines) == 39
    assert len(isa.nodes) == 1729
    assert len(isa.nodechains) == 98
    assert len(isa.elements) == 822
    assert len(isa.superelements) == 3
    assert len(isa.get_subregion('asdf').elements()) == 756


def test_no_such_subregion(isa):
    with pytest.raises(ValueError) as excinfo:
        n = 'foo'
        isa.get_subregion(n)
    assert str(excinfo.value) == (
        'no such subregion "{}" in this model'.format(n))


def test_points(isa):
    p = isa.points[0].x, isa.points[0].y
    p_xy = isa.points[0].xy

    assert p == p_xy == pytest.approx((0.008339, 0.011702), abs=1e-5)


def test_lines(isa):
    ln = isa.lines[0]

    assert type(ln.p1) == isa7.Point
    assert type(ln.p2) == isa7.Point
    assert ln.p1.xy == pytest.approx((0.008339, 0.011702), abs=1e-5)
    assert ln.p2.xy == pytest.approx((0.007388, 0.007190), abs=1e-5)


def test_nodes(isa):
    nd = isa.nodes[-1]

    assert type(nd) == isa7.Node
    assert nd.xy == (nd.x, nd.y) == pytest.approx(
        (0.008787, 0.007736), abs=1e-5)


def test_nodechains(isa):
    nc = isa.nodechains[0]

    assert type(nc) == isa7.NodeChain
    assert type(nc.nodes[0]) == isa7.Node


def test_nodechains_reverse(isa):
    for nc in isa.nodechains:
        assert nc.reverse().nodes == nc.nodes[::-1]


def test_elements(isa):
    el = isa.elements[0]

    assert type(el) == isa7.Element
    for v in el.vertices:
        assert type(v) == isa7.Node
    assert el.vertices[0].xy == pytest.approx((0.002202, 0.003197),
                                              abs=1e-5)


def test_superelements(isa):
    se = isa.superelements[0]

    assert type(se) == isa7.SuperElement
    for el in se.elements:
        assert type(el) == isa7.Element
    for nc in se.nodechains:
        assert type(nc) == isa7.NodeChain
        assert nc.key == se.nc_keys[se.nodechains.index(nc)]

        
@pytest.fixture
def disp_stat():
    filename = 'tests/data/test_disp_stat.ISA7'
    return isa7.read(filename)


def test_subregions(disp_stat):
    sr = disp_stat.subregions[0]

    assert type(sr) == isa7.SubRegion
    assert sr.name == "Stat"
    for se in sr.superelements:
        assert type(se) == isa7.SuperElement
    for nc in sr.nodechains:
        assert type(nc) == isa7.NodeChain


def test_windings(disp_stat):
    wd = disp_stat.windings[0]
    
    assert type(wd) == isa7.Winding
    assert wd.name == "Stra"
    for sr in wd.subregions:
        assert type(sr) == isa7.SubRegion
    assert wd.num_turns == 100
