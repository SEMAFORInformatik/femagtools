import pytest
from femagtools import vtu


@pytest.fixture
def vtu_data():
    filename = 'tests/data/zzz_pm_model_ts_results_1/zzz_pm_model_ts_0000.vtu'
    return vtu.read(filename)


def test_read(vtu_data):
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
