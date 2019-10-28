import meshio
from femagtools import convert, isa7


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

    m = meshio.read(msh)
    assert len(m.points) == 2340
    assert len(m.point_data["potential"]) == len(m.points)
    assert len(m.cells["triangle"]) == 1522
    assert len(m.cells["quad"]) == 1506
    assert len(m.field_data.keys()) == 28
    num_cells = sum(map(len, m.cells.values()))
    for data in m.cell_data["triangle"]:
        len_data = sum([len(m.cell_data[ctype][data])
                        for ctype in m.cell_data])
        assert len_data == num_cells


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
    assert len(mesh.cells["triangle"]) == 2


def test_nastran_quads(tmpdir):
    msh = str(tmpdir.join("quads.msh"))

    convert.to_msh("tests/data/quads.nas", msh)

    mesh = meshio.read(msh)

    assert len(mesh.points) == 6
    assert len(mesh.cells["quad"]) == 2


def test_nastran_superelements(tmpdir):
    msh = str(tmpdir.join("superelements.msh"))

    convert.to_msh("tests/data/superelements.nas", msh)

    mesh = meshio.read(msh)

    assert len(mesh.points) == 6
    assert len(mesh.cells["triangle"]) == 2
    assert len(mesh.cells["quad"]) == 1

    assert mesh.cell_data["triangle"]["gmsh:geometrical"][0] == 2
    assert mesh.cell_data["triangle"]["gmsh:geometrical"][1] == 1

    assert mesh.cell_data["quad"]["gmsh:geometrical"][0] == 2


def test_vtu(tmpdir):
    vtu = str(tmpdir.join("magnsec.vtu"))
    vtu2 = str(tmpdir.join("magnsec2.vtu"))

    convert.to_vtu("tests/data/magnsec.ISA7", vtu)

    isa = isa7.read("tests/data/magnsec.ISA7")
    convert.to_vtu(isa, vtu2)

    with open(vtu) as f:
        file1 = f.read()
    with open(vtu2) as f:
        file2 = f.read()
    assert file1 == file2
        
    m = meshio.read(vtu)
    assert len(m.points) == 2340
    assert len(m.point_data["potential"]) == len(m.points)
    assert len(m.cells["triangle"]) == 1522
    assert len(m.cells["quad"]) == 1506
    num_cells = sum(map(len, m.cells.values()))
    for data in m.cell_data["triangle"]:
        len_data = sum([len(m.cell_data[ctype][data])
                        for ctype in m.cell_data])
        assert len_data == num_cells


def test_vtu_triangles(tmpdir):
    vtu = str(tmpdir.join("triangles.vtu"))
    convert.to_vtu("tests/data/triangles.ISA7", vtu)
    with open(vtu) as f:
        assert len(f.readlines()) == 28


def test_vtu_quads(tmpdir):
    vtu = str(tmpdir.join("quads.vtu"))
    convert.to_vtu("tests/data/quads.ISA7", vtu)
    with open(vtu) as f:
        assert len(f.readlines()) == 28


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
