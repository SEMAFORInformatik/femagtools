import femagtools.asm
import pathlib
import pytest

expected = {'version': '9.2.x-202-geff049a1',
            'project': './LS4p.nc', 'filename': 'TEST_001', 'date': '2021-05-27T13:23',
            'wdgconn': 'star',
            'f1': [50.0, 50.0, 50.0], 'r1': 1.2, 'lfe': 160.0,
            'rbarlen': 246.0, 'num_phases': 3, 'p': 2, 'p_gen': 1,
            'num_par_wdgs': 1, 'mcfile': 'M330-50A',
            'felosscoeff': 1.44, 'maxiters': 300,
            'permchg': 0.0, 's': [0.0, 0.025, 0.05], 'T': [0.0, 48.8, 92.4],
            'u1': [380.0, 380.0, 380.0],
            'Tp2': [0.0, 45.5, 84.7], 'p2': [0.0, 179.0, 665.0], 'pfe1': [],
            'i1': [8.52, 14.9, 25.3], 'p1': [0.0, 7230.0, 13600.0],
            'cosphi': [0.0, 0.738, 0.818],
            'pcu': [261.32544, 799.236, 2304.324],
            'pltotal': [261.32544, 978.236, 2969.324],
            'lh': 0.079515, 'ls1': 0.00245,
            'ls2': 0.00925, 'r2': 0.406538, 'sk': 0.1135}


def test_read_asm():
    datadir = pathlib.Path(__file__).resolve().parent.joinpath('data')

    content = (datadir / 'test.ASM').read_text()
    r = femagtools.asm.read(content.split('\n'))
    assert r.keys() == expected.keys()
    for k in r.keys():
        if isinstance(r[k], float) or isinstance(r[k], list):
            assert pytest.approx(r[k], rel=0.01) == expected[k]
        else:
            assert r[k] == expected[k]
