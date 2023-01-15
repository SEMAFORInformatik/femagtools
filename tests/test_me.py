from femagtools import me

def test_me():
    r = me.get_eigenvectors(['tests/data/Mode_001.txt'], [])
    assert r[2]['0']['eigenvecs'].shape[0] == 5796
    assert round(r[1][0], 2) == 52.07
