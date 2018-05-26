
from femagtools import forcedens


def test_plt_read():
    filename = 'tests/data/PLT.0'
    fdens = forcedens.ForceDensity()
    fdens.read(filename)

    assert fdens.title == 'Load PM-Syn_motor airgap 1, No Skewing'
    assert len(fdens.positions) == 31
    assert sorted(fdens.positions[0].keys()) == ['B_N', 'B_T', 'FN', 'FT',
                                                 'Radius', 'X',
                                                 'column_units', 'position',
                                                 'unit']

