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
    assert {'r2': 0.00011541,
            'ls2': 3.57991e-07,
            'zeta2': 0.9891359043215603,
            'pl2v': 0.9130062527810479} == r
