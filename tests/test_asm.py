import femagtools.asm
import pathlib

expected = {'version': '9.2.x-202-geff049a1', 'project': './LS4p.nc', 'filename': 'TEST_001', 'date': '2021-05-27T13:23',
            'u1': 219.0, 'wdgconn': 'star', 'f1': [50.0, 50.0, 50.0], 'r1': 1.2, 'xs1': 0.77, 'lfe': 160.0,
            'rbarlen': 246.0, 'num_phases': 3, 'p': 2, 'p_gen': 1, 'num_par_wdgs': 1, 'mcfile': 'M330-50A',
            'felosscoeff': '1.44', 'maxiters': 300, 'permchg': 0.0, 's': [0.0, 2.5, 5.0], 'T': [0.0, 48.8, 92.4],
            'un': [380.0, 380.0, 380.0], 'i1': [8.52, 14.9, 25.3], 'p1': [0.0, 7.23, 13.6], 'cosphi': [0.0, 0.738, 0.818]}


def test_read_asm():
    datadir = pathlib.Path(__file__).resolve().parent.joinpath('data')

    content = (datadir / 'test.ASM').read_text()
    r = femagtools.asm.read(content)
    assert r == expected
