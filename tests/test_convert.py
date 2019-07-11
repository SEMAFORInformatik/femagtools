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



