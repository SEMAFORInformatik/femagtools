#!/usr/bin/env python
#
import pytest
import femagtools.windings

wpar = dict(Q = 54, p = 6, m = 3,
            windings = {1: {
                'dir': [1, 1, 1, -1, -1, -1],
                'N': [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
                'PHI': [3.3333, 3.3333, 10.0, 30.0, 36.6666, 36.6666]}})


@pytest.fixture()
def wdg():
    return femagtools.windings.Windings(wpar)


def test_slots(wdg):
    assert wdg.slots(1).tolist() == [
        [ 0,  0,  1,  4,  5,  5],
        [ 9,  9, 10, 13, 14, 14],
        [18, 18, 19, 22, 23, 23],
        [27, 27, 28, 31, 32, 32],
        [36, 36, 37, 40, 41, 41],
        [45, 45, 46, 49, 50, 50]]

def test_axis(wdg):
    assert round(wdg.axis(), 3) == 0.349
