import femagtools.machine
from femagtools.job import Task
import pathlib
import pytest


@pytest.fixture
def data_dir():
    return pathlib.Path(__file__).with_name('data')


def test_ecsim_eval(data_dir):
    task = Task(id=1, directory=data_dir / 'parident')
    eval = femagtools.machine.im._eval_ecsim()
    r = eval(task)
    expected = {'r2': 0.00011541,
                'ls2': 3.57991e-07,
                'zeta2': 0.9891,
                'pl2v': 0.913}
    assert expected['r2'] == r['r2']
    assert expected['ls2'] == r['ls2']
    assert expected['zeta2'] == round(r['zeta2'], 4)
    assert expected['pl2v'] == round(r['pl2v'], 3)
