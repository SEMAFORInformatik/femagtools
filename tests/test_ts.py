import pytest
from femagtools import ts


@pytest.fixture
def losses():
    dirname = 'tests/data/zzz_pm_model_ts_results_1'
    modelname = 'tests/data/zzz_pm_model_ts'
    return ts.Losses(modelname, dirname)


def test_read(losses):
    assert [{'key': 1, 'name': 'stfe', 'losses': 0.0},
            {'key': 2, 'name': 'rofe', 'losses': 0.0},
            {'key': 3, 'name': 'wefe', 'losses': 0.0},
            {'key': 4, 'name': '', 'losses': 0.001912685411090982},
            {'key': 5, 'name': '', 'losses': 0.0019698298472865085},
            {'key': 6, 'name': '', 'losses': 0.0018862924108952414}] == losses.ohm_lossenergy(0.0, 0.0)
    assert [{'key': 1, 'name': 'stfe', 'losses': 0.0},
            {'key': 2, 'name': 'rofe', 'losses': 0.0},
            {'key': 3, 'name': 'wefe', 'losses': 0.0},
            {'key': 4, 'name': '', 'losses': 0.09808646181380988},
            {'key': 5, 'name': '', 'losses': 0.10101694663178479},
            {'key': 6, 'name': '', 'losses': 0.09673297420375133}] == losses.ohm_powerlosses()
    assert [{'key': 1, 'name': 'stfe', 'losses': 0.0},
            {'key': 2, 'name': 'rofe', 'losses': 0.0},
            {'key': 3, 'name': 'wefe', 'losses': 0.0},
            {'key': 4, 'name': '', 'losses': 0.10108592471259238},
            {'key': 5, 'name': '', 'losses': 0.10108593475500291},
            {'key': 6, 'name': '', 'losses': 0.10108575785653288}] == losses.ohm_powerlosses_fft()
