import pytest
from femagtools import nc
from femagtools import isa7


@pytest.fixture
def model():
    filename = 'src/tests/data/minimal.nc'
    return nc.read(filename)

@pytest.fixture
def pm():
    filename = 'src/tests/data/pm_data.nc'
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

def test_airgap_center_elements(pm):
    assert len(pm.airgap_center_elements) == 180

def test_areas(pm):
    a = pm.get_areas()
    assert a[0]['iron'] == pytest.approx(0.003550, abs=1e-5)
    assert a[0]['slots'] == pytest.approx(0.00184, abs=1e-5)
    assert a[0]['magnets'] == 0.0
    assert a[1]['iron'] == pytest.approx(0.0007222, abs=1e-5)
    assert a[1]['slots'] == 0
    assert a[1]['magnets'] == pytest.approx(0.000345, abs=1e-5)

def test_geom(pm):
    mag_spels = pm.magnet_super_elements()
    assert len(mag_spels) == 5
    g = mag_spels[0].get_rect_geom()
    assert g['w'] == pytest.approx(0.0112, abs=1e-4)
    assert g['h'] == pytest.approx(0.00308, abs=1e-5)
    assert g['cxy'] == pytest.approx((0.02317, 0.007528), abs=1e-5)
    assert g['area'] == pytest.approx(3.45e-05, abs=1e-6)
    assert g['alpha'] == pytest.approx(1.885, abs=1e-3)

def test_calc_iron_loss(pm):
    import numpy as np
    def pfe(Bxnu, Bynu, fnu, losscoeffs):
        basfrq = 50
        basind = 1.5
        b21 = np.linalg.norm((Bxnu, Bynu), axis=0)
        b2 = (b21/basind)
        f = fnu/basfrq
        hch = f
        hcw = f**2
        phy = f*b2
        pec = f**2*b2**2
        pex = f**1.5*b2**1.5
        return phy, pec, pex

    icur = 0
    ibeta = 0
    pfe = pm.calc_iron_loss(icur, ibeta, pfe)
    assert pfe['Rüc'] == [0, 0, 0]
    assert pfe['Stat'] == pytest.approx([5.61, 46.982, 15.92], abs=1e-2)


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
    filename = 'src/tests/data/test_disp_stat.nc'
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
