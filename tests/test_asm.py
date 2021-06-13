import femagtools.asm
import pathlib

expected = {'version': '9.2.x-202-geff049a1',
            'project': './LS4p.nc', 'filename': 'TEST_001', 'date': '2021-05-27T13:23',
            'u1': 219.0, 'wdgconn': 'star',
            'f1': [50.0, 50.0, 50.0], 'r1': 1.2, 'xs1': 0.77, 'lfe': 160.0,
            'rbarlen': 246.0, 'num_phases': 3, 'p': 2, 'p_gen': 1,
            'num_par_wdgs': 1, 'mcfile': 'M330-50A',
            'felosscoeff': '1.44', 'maxiters': 300,
            'permchg': 0.0, 's': [0.0, 2.5, 5.0], 'T': [0.0, 48.8, 92.4],
            'u1': [380.0, 380.0, 380.0],
            'Tp2': [0.0, 45.5, 84.7], 'p2': [0.0, 179.0, 665.0], 'pfe1': [],
            'i1': [8.52, 14.9, 25.3], 'p1': [0.0, 7230.0, 13600.0],
            'cosphi': [0.0, 0.738, 0.818],
            'pcu': [261.32543999999996, 799.2360000000001, 2304.324],
            'pltotal': [261.32543999999996, 978.2360000000001, 2969.324]}


def test_read_asm():
    datadir = pathlib.Path(__file__).resolve().parent.joinpath('data')

    content = (datadir / 'test.ASM').read_text()
    r = femagtools.asm.read(content.split('\n'))
    assert r == expected
