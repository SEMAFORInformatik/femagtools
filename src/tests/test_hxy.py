import femagtools.hxy
import pathlib
import pytest


@pytest.fixture
def data_dir():
    return pathlib.Path(__file__).with_name('data') / 'hxy'


def test_read(data_dir):
    num_magnets = 2
    magnets = femagtools.hxy.read(data_dir / 'PM270L8_011.hxy',
                                  num_magnets)
    assert len(magnets) == num_magnets
    assert set([len(m['e'][0]) for m in magnets]) == {169, 189}
    assert [m['pos'][0] for m in magnets] == [0.0, 0.0]
    assert sorted([m['havg'][0] for m in magnets]) == pytest.approx([156.3, 156.5], 0.1)
    assert sorted([m['hmax'][0] for m in magnets]) == pytest.approx([195.6, 304.4], 0.1)
