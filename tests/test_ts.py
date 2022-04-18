import pytest
from femagtools import ts


@pytest.fixture
def losses():
    dirname = 'tests/data/zzz_pm_model_ts_results_1'
    modelname = 'tests/data/zzz_pm_model_ts'
    return ts.Losses(modelname, dirname)


def test_read(losses):
    e = losses.ohm_lossenergy(0.0, 0.0)
    assert len(e) == 6
    assert {'key': 1, 'name': 'stfe', 'losses': 0.0} == e[0]
    assert {'key': 2, 'name': 'rofe', 'losses': 0.0} == e[1]
    assert {'key': 3, 'name': 'wefe', 'losses': 0.0} == e[2]
    assert pytest.approx(
        {'key': 4, 'name': '', 'losses': 0.001913}, abs=1e-5) == e[3]
    assert pytest.approx(
        {'key': 5, 'name': '', 'losses': 0.001967}, abs=1e-5) == e[4]
    assert pytest.approx(
        {'key': 6, 'name': '', 'losses': 0.001886}, abs=1e-5) == e[5]
    p = losses.ohm_powerlosses()
    assert len(p) == 6
    assert {'key': 1, 'name': 'stfe', 'losses': 0.0} == p[0]
    assert {'key': 2, 'name': 'rofe', 'losses': 0.0} == p[1]
    assert {'key': 3, 'name': 'wefe', 'losses': 0.0} == p[2]
    assert pytest.approx(
        {'key': 4, 'name': '', 'losses': 0.09809}, abs=1e-5) == p[3]
    assert pytest.approx(
        {'key': 5, 'name': '', 'losses': 0.10102}, abs=1e-5) == p[4]
    assert pytest.approx(
        {'key': 6, 'name': '', 'losses': 0.09673}, abs=1e-5) == p[5]
    p = losses.ohm_powerlosses_fft()
    assert len(p) == 6
    assert {'key': 1, 'name': 'stfe', 'losses': 0.0} == p[0]
    assert {'key': 2, 'name': 'rofe', 'losses': 0.0} == p[1]
    assert {'key': 3, 'name': 'wefe', 'losses': 0.0} == p[2]
    assert pytest.approx(
        {'key': 4, 'name': '', 'losses': 0.101086}, abs=1e-5) == p[3]
    assert pytest.approx(
        {'key': 5, 'name': '', 'losses': 0.101086}, abs=1e-5) == p[4]
    assert pytest.approx(
        {'key': 6, 'name': '', 'losses': 0.101086}, abs=1e-5) == p[5]
