import pathlib
from femagtools.dxfsl.converter import convert

def test_convert():
    p = pathlib.Path(__file__).parent / 'data' / 'IPM-130-4.dxf'
    r = convert(str(p))

    assert r['num_poles'] == 4
    assert r['tot_num_slot'] == 12
    totnumsl = [l for l in r['fsl'] if l.startswith('m.tot_num_slot')]
    assert len(totnumsl) == 1
    assert totnumsl[0].split()[-1] == '12'
