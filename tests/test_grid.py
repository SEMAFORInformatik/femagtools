#!/usr/bin/env python
#
import femagtools.grid


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


def test_baskets():
    x = list(range(5))*5
    baskets = femagtools.grid.baskets(x, 5)

    for i, p in enumerate(baskets):
        assert p == [0, 1, 2, 3, 4]
    assert i == 4

        
