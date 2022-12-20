import numpy as np
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


def test_noloadrot_eval(data_dir):
    task = Task(id=1, directory=data_dir / 'parident')
    i1tab = [2.0384197650236815, 4.076839530047363,
             6.115259295071045, 8.153679060094726, 10.192098825118407]
    u1ph = 400/np.sqrt(3)
    eval = femagtools.machine.im._eval_noloadrot()
    r = eval(task)
    i0 = [2.0384,  4.0768,  6.1153,  8.1537, 10.1921]
    psi1_0 = [0.2725, 0.5301, 0.7106, 0.8091, 0.8641]
    Bamp = [0.3091, 0.6013, 0.8058, 0.9159, 0.9769]
    assert i0 == [round(x, 4) for x in r['i1_0']]
    assert psi1_0 == [round(x, 4) for x in np.mean(r['psi1_0'], axis=1)]
    assert Bamp == [round(x, 4) for x in np.mean(r['Bamp'], axis=1)]
