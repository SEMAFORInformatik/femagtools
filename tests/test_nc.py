import pytest
from femagtools import nc
from femagtools import isa7


@pytest.fixture
def model():
    filename = 'tests/data/minimal.nc'
    return nc.read(filename)


def test_read(model):
    assert len(model.points) == 39
    assert len(model.lines) == 39 
    assert model.points[0].x == pytest.approx(0.008339,
                                            abs=1e-5)
    

def test_objects(model):
    assert len(model.nodes) == 1729
    assert len(model.nodechains) == 98
    assert len(model.elements) == 822
    assert len(model.superelements) == 3
    assert len(model.get_subregion('asdf').elements()) == 756

def test_no_such_subregion(model):
    with pytest.raises(ValueError) as excinfo:
        n = 'foo'
        model.get_subregion(n)
    assert str(excinfo.value) == (
        'no such subregion "{}" in this model'.format(n))

def test_points(model):
    p = model.points[0].x, model.points[0].y
    p_xy = model.points[0].xy

    assert p == p_xy == pytest.approx((0.008339, 0.011702), abs=1e-5)


def test_lines(model):
    ln = model.lines[0]

    assert type(ln.p1) == isa7.Point
    assert type(ln.p2) == isa7.Point
    assert ln.p1.xy == pytest.approx((0.008339, 0.011702), abs=1e-5)
    assert ln.p2.xy == pytest.approx((0.007388, 0.007190), abs=1e-5)


def test_nodes(model):
    nd = model.nodes[-1]

    assert type(nd) == isa7.Node
    assert nd.xy == (nd.x, nd.y) == pytest.approx(
        (0.008787, 0.007736), abs=1e-5)


def test_nodechains(model):
    nc = model.nodechains[0]

    assert type(nc) == isa7.NodeChain
    assert type(nc.nodes[0]) == isa7.Node


def test_nodechains_reverse(model):
    for ndc in model.nodechains:
        assert ndc.reverse().nodes == ndc.nodes[::-1]


def test_elements(model):
    el = model.elements[0]

    assert type(el) == isa7.Element
    for v in el.vertices:
        assert type(v) == isa7.Node
    assert el.vertices[0].xy == pytest.approx((0.002202, 0.003197),
                                              abs=1e-5)


def test_superelements(model):
    se = model.superelements[0]

    assert type(se) == isa7.SuperElement
    for el in se.elements:
        assert type(el) == isa7.Element
    for ndc in se.nodechains:
        assert type(ndc) == isa7.NodeChain
        assert ndc.key == se.nc_keys[se.nodechains.index(ndc)]

        
@pytest.fixture
def disp_stat():
    filename = 'tests/data/test_disp_stat.nc'
    return nc.read(filename)


def test_subregions(disp_stat):
    sr = disp_stat.subregions[0]

    assert type(sr) == isa7.SubRegion
    assert sr.name == "Stat"
    for se in sr.superelements:
        assert type(se) == isa7.SuperElement
    for ndc in sr.nodechains:
        assert type(ndc) == isa7.NodeChain


def test_windings(disp_stat):
    wd = disp_stat.windings[0]
    
    assert type(wd) == isa7.Winding
    assert wd.name == "Stra"
    for sr in wd.subregions:
        assert type(sr) == isa7.SubRegion
    assert wd.num_turns == 100


def test_magnet_super_elements(disp_stat):
    sekeys = [se.key for se in disp_stat.magnet_super_elements()]
    assert sekeys == [98, 101, 74, 77, 80, 83, 86, 89, 92, 95]
