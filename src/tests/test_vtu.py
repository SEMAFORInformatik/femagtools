import pytest
import pathlib
from femagtools import vtu


@pytest.fixture
def ts_data_dir():
    return pathlib.Path(__file__).with_name('data') / 'zzz_pm_model_ts_results_1/'

@pytest.fixture
def demag_data_dir():
    return pathlib.Path(__file__).with_name('data') / 'demag-vtu/'


def test_read(ts_data_dir):
    vtu_data = vtu.read(ts_data_dir / 'zzz_pm_model_ts_0000.vtu')
    assert vtu_data.field_data_names == [
        'time [s]', 'angle [rad]',
        'speed [rad/s]', 'torque [Nm]',
        'winding 1 current [A]', 'winding 1 voltage [V]',
        'winding 2 current [A]', 'winding 2 voltage [V]',
        'winding 3 current [A]', 'winding 3 voltage [V]', 'branch names']
    assert vtu_data.point_data_names == ['vector potential']
    assert vtu_data.cell_data_names == [
        'b', 'curd', 'demagnetization', 'current [A]', 'voltage [V]']

    b = vtu_data.get_data_vector('b', 1)
    assert b[0] == pytest.approx([-0.003114], abs=1e-5)
    assert b[1] == pytest.approx([-0.00313], abs=1e-5)
    assert b[2] == [0.0]

def test_demag(demag_data_dir):
    vtu_data = vtu.read(demag_data_dir / 'PM_130_L10_0000.vtu')
    keys = [7412, 7413, 7414, 7415, 7416]
    class Element:
        def __init__(self, key):
            self.key = key
    elements = [Element(k) for k in keys]
    expected = [241.5, 250.2, 262.3, 276.5, 390.5]
    actual = vtu_data.demag(elements)
    assert len(actual) == 1
    assert actual[0] == pytest.approx(expected, abs=0.1)
