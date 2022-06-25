import meshio
from femagtools import convert, isa7
import xml.etree.ElementTree as ET


def meshio_major_version():
    return int(meshio.__version__.split('.')[0])


def test_msh(tmpdir):
    msh = str(tmpdir.join("magnsec.msh"))
    msh2 = str(tmpdir.join("magnsec2.msh"))

    convert.to_msh("tests/data/magnsec.ISA7", msh)

    isa = isa7.read("tests/data/magnsec.ISA7")
    convert.to_msh(isa, msh2)

    with open(msh) as f:
        file1 = f.read()
    with open(msh2) as f:
        file2 = f.read()
    assert file1 == file2
    # meshio 5 bug https://github.com/nschloe/meshio/issues/1257
    if meshio_major_version() < 5:
        m = meshio.read(msh)
        assert len(m.points) == 2340
        assert len(m.point_data["potential"]) == len(m.points)
        assert len(m.cells[1][1]) == 1522
        assert len(m.cells[2][1]) == 1506
        assert len(m.field_data.keys()) == 28


def test_msh_triangles(tmpdir):
    msh = str(tmpdir.join("triangles.msh"))
    convert.to_msh("tests/data/triangles.ISA7", msh)


def test_msh_quads(tmpdir):
    msh = str(tmpdir.join("quads.msh"))
    convert.to_msh("tests/data/quads.ISA7", msh)


def test_nastran_float():

    for seven in ["7.0", ".7E1", "0.7+1", ".70+1", "7.E+0", "70.-1"]:
        assert convert._nastran_real_to_float(seven) == 7.0

    for negative_seven in ["-7.0", "-.7E1", "-0.7+1", "-.70+1", "-7.E+0", "-70.-1"]:
        assert convert._nastran_real_to_float(negative_seven) == -7.0


def test_nastran_triangles(tmpdir):
    msh = str(tmpdir.join("triangles.msh"))

    convert.to_msh("tests/data/triangles.nas", msh)

    mesh = meshio.read(msh)

    assert len(mesh.points) == 4
    if meshio_major_version() < 5:
        assert mesh.cells[0][0] == "triangle" and len(mesh.cells[0][1]) == 2
    else:
        assert mesh.cells[0].type == "triangle" and mesh.cells[0].data.shape == (
            2, 3)


def test_nastran_quads(tmpdir):
    msh = str(tmpdir.join("quads.msh"))

    convert.to_msh("tests/data/quads.nas", msh)

    mesh = meshio.read(msh)

    assert len(mesh.points) == 6
    if meshio_major_version() < 5:
        assert mesh.cells[0][0] == "quad" and len(mesh.cells[0][1]) == 2
    else:
        assert mesh.cells[0].type == "quad" and mesh.cells[0].data.shape == (
            2, 4)


def test_nastran_superelements(tmpdir):
    msh = str(tmpdir.join("superelements.msh"))

    convert.to_msh("tests/data/superelements.nas", msh)

    mesh = meshio.read(msh)

    assert len(mesh.points) == 6
    if meshio_major_version() < 5:
        assert mesh.cells[0][0] == "triangle" and len(mesh.cells[0][1]) == 2
        assert mesh.cells[1][0] == "quad" and len(mesh.cells[1][1]) == 1
    else:
        assert mesh.cells[0].type == "triangle"
        assert mesh.cells[0].data.shape == (2, 3)
        assert mesh.cells[1].type == "quad"
        assert mesh.cells[1].data.shape == (1, 4)

    assert mesh.cell_data["gmsh:geometrical"][0][0] == 2
    assert mesh.cell_data["gmsh:geometrical"][0][1] == 1

    assert mesh.cell_data["gmsh:geometrical"][1][0] == 2


def test_vtu(tmpdir):
    vtu = str(tmpdir.join("magnsec.vtu"))
    convert.to_vtu("tests/data/magnsec.ISA7", vtu)

    tree = ET.parse(vtu)
    root = tree.getroot()
    #
    assert [child.tag for child in root] == ['UnstructuredGrid']
    assert [child.tag for child in root[0]] == ['Piece']
    assert [child.tag for child in root[0][0]] == ['Points',
                                                   'Cells',
                                                   'PointData',
                                                   'CellData']
    assert [child.tag for child in root[0][0][1]] == ['DataArray',
                                                      'DataArray',
                                                      'DataArray']

    m = meshio.read(vtu)
    assert len(m.points) == 2340
    assert len(m.point_data["potential"]) == len(m.points)
    if meshio_major_version() < 5:
        assert m.cells[1][0] == "triangle" and len(m.cells[1][1]) == 1522
        assert m.cells[2][0] == "quad" and len(m.cells[2][1]) == 1506
    else:
        assert m.cells[1].type == "triangle"
        assert m.cells[1].data.shape == (1522, 3)
        assert m.cells[2].type == "quad"
        assert m.cells[2].data.shape == (1506, 4)


def test_vtu_triangles(tmpdir):
    vtu = str(tmpdir.join("triangles.vtu"))
    convert.to_vtu("tests/data/triangles.ISA7", vtu)
    tree = ET.parse(vtu)
    root = tree.getroot()
    assert [child.tag for child in root] == ['UnstructuredGrid']
    assert [child.tag for child in root[0]] == ['Piece']
    assert [child.tag for child in root[0][0]] == ['Points',
                                                   'Cells',
                                                   'PointData',
                                                   'CellData']
    assert [child.tag for child in root[0][0][1]] == ['DataArray',
                                                      'DataArray',
                                                      'DataArray']


def test_vtu_quads(tmpdir):
    vtu = str(tmpdir.join("quads.vtu"))
    convert.to_vtu("tests/data/quads.ISA7", vtu)
    tree = ET.parse(vtu)
    root = tree.getroot()
    assert [child.tag for child in root] == ['UnstructuredGrid']
    assert [child.tag for child in root[0]] == ['Piece']
    assert [child.tag for child in root[0][0]] == ['Points',
                                                   'Cells',
                                                   'PointData',
                                                   'CellData']
    assert [child.tag for child in root[0][0][1]] == ['DataArray',
                                                      'DataArray',
                                                      'DataArray']


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
