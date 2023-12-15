import femagtools.hxy
import pathlib
import pytest


@pytest.fixture
def data_dir():
    return pathlib.Path(__file__).with_name('data')


def test_read(data_dir):
    num_magnets = 2
    magnets = femagtools.hxy.read(str(data_dir / 'PM270L8_011.hxy'), num_magnets)
    assert len(magnets) == num_magnets
