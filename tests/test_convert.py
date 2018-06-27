import pytest
import meshio
from femagtools import convert, isa7


def test_msh(tmpdir):
    msh = str(tmpdir.join("magnsec.msh"))
    msh2 = str(tmpdir.join("magnsec2.msh"))

    convert.to_msh("tests/data/magnsec.ISA7", msh)

    isa = isa7.read("tests/data/magnsec.ISA7")
    convert.to_msh(isa, msh2)

    with open(msh) as f:
        assert len(f.readlines()) == 8082
    with open(msh2) as f:
        assert len(f.readlines()) == 8082

    mesh = meshio.read(msh)
    mesh2 = meshio.read(msh2)
    for m in mesh, mesh2:
        assert len(m.points) == 2340
        assert len(m.point_data["potential"]) == 2340
        assert len(m.cells["triangle"]) == 1522
        assert len(m.cells["quad"]) == 1506
        assert len(m.field_data.keys()) == 28


def test_vtu(tmpdir):
    vtu = str(tmpdir.join("magnsec.vtu"))
    vtu2 = str(tmpdir.join("magnsec2.vtu"))

    convert.to_vtu("tests/data/magnsec.ISA7", vtu)

    isa = isa7.read("tests/data/magnsec.ISA7")
    convert.to_vtu(isa, vtu2)

    with open(vtu) as f:
        assert len(f.readlines()) == 22
    with open(vtu2) as f:
        assert len(f.readlines()) == 22

    mesh = meshio.read(vtu)
    mesh2 = meshio.read(vtu2)
    for m in mesh, mesh2:
        assert len(m.points) == 2340
        assert len(m.point_data["potential"]) == 2340
        assert len(m.cells["triangle"]) == 1522
        assert len(m.cells["quad"]) == 1506


def test_geo(tmpdir):
    geo = str(tmpdir.join("magnsec.geo"))
    convert.to_geo("tests/data/magnsec.ISA7", geo)
    with open(geo) as f:
        assert len(f.readlines()) == 2537


def test_geo_extrude(tmpdir):
    geo = str(tmpdir.join("magnsec.geo"))
    convert.to_geo("tests/data/magnsec.ISA7", geo, 0.01, 3, True)
    with open(geo) as f:
        assert len(f.readlines()) == 2981



