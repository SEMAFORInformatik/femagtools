#!/usr/bin/env python
#
import pytest
import femagtools.grid
import femagtools.parstudy
import numpy as np
import functools


def test_create_parameter_range():
    x = [(1, 2, 3), (4, 5), (6, 7)]

    r = femagtools.grid.create_parameter_range(x)

    assert r.tolist() == [[1, 4, 6],
                          [2, 4, 6],
                          [3, 4, 6],
                          [1, 5, 6],
                          [2, 5, 6],
                          [3, 5, 6],
                          [1, 4, 7],
                          [2, 4, 7],
                          [3, 4, 7],
                          [1, 5, 7],
                          [2, 5, 7],
                          [3, 5, 7]]


@pytest.fixture
def decision_vars():
    return [
        {"steps": 3, "bounds": [-50, 0],
         "name": "angl_i_up", "label":"Beta"},
        {"steps": 3, "bounds": [100, 200],
         "name": "current", "label":"Current/A"}
    ]


def test_baskets():
    x = list(range(5))*5
    baskets = femagtools.parstudy.baskets(x, 5)

    for i, p in enumerate(baskets):
        assert p == [0, 1, 2, 3, 4]
    assert i == 4


def test_report(decision_vars):
    objective_vars = [
        {"name": "dqPar.torque[-1]",
         "label": "Load Torque/Nm"},
        {"name": "torque[-1].ripple",
         "label": "Torque Ripple/Nm"},
        {"name": "machine.plfe[-1]",
         "label": "Iron Losses/W"}
    ]

    shape = [len(objective_vars)] + [d['steps'] for d in decision_vars]
    a = [0]*len(objective_vars)*functools.reduce(
        (lambda x, y: x * y), [d['steps'] for d in decision_vars])
    objectives = np.reshape(a, shape)

    expected = [[d['label'] for d in decision_vars] +
                [o['label'] for o in objective_vars] + ['Directory'],
                [d['name'] for d in decision_vars] +
                [o['name'] for o in objective_vars],
                [-50.0, 100.0, 0.0, 0.0, 0.0, 0],
                [-25.0, 100.0, 0.0, 0.0, 0.0, 1],
                [-0.0,  100.0, 0.0, 0.0, 0.0, 2],
                [-50.0, 150.0, 0.0, 0.0, 0.0, 3],
                [-25.0, 150.0, 0.0, 0.0, 0.0, 4],
                [-0.0,  150.0, 0.0, 0.0, 0.0, 5],
                [-50.0, 200.0, 0.0, 0.0, 0.0, 6],
                [-25.0, 200.0, 0.0, 0.0, 0.0, 7],
                [-0.0, 200.0, 0.0, 0.0, 0.0, 8]]

    domain = [list(np.linspace(d['bounds'][0], d['bounds'][1], d['steps']))
              for d in decision_vars]
    par_range = femagtools.grid.create_parameter_range(domain)
    assert expected == femagtools.parstudy.get_report(decision_vars,
                                                      objective_vars, objectives,
                                                      par_range)


def test_sobol(decision_vars):
    parvar = femagtools.parstudy.Sobol('.')
    N = 4
    n, d, r = parvar._get_names_and_range(decision_vars, N)
    assert [d['name'] for d in decision_vars] == n
    assert r.shape == (N, len(decision_vars))


def test_lhs(decision_vars):
    parvar = femagtools.parstudy.LatinHypercube('.')
    N = 4
    n, d, r = parvar._get_names_and_range(decision_vars, N)
    assert [d['name'] for d in decision_vars] == n
    assert r.shape == (N, len(decision_vars))
