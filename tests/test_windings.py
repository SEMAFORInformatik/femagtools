#!/usr/bin/env python
#
import pytest
import femagtools.windings


@ pytest.fixture()
def wdg():
    wpar = dict(Q=54, p=6, m=3,
                windings={
                    1: {'dir': [1, 1, 1, -1, -1, -1],
                        'N': [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
                        'R': [1, 0, 1, 0, 1, 0],
                        'PHI': [3.3333, 3.3333, 10.0, 30.0, 36.6666, 36.6666]},
                    2: {'dir': [1, 1, 1, -1, -1, -1],
                        'N': [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
                        'R': [1, 0, 1, 0, 1, 0],
                        'PHI': [23.333333333333332,
                                23.333333333333332,
                                30.0,
                                50.00000000000001,
                                56.66666666666667,
                                56.66666666666667]},
                    3: {'dir': [-1, -1, -1, 1, 1, 1],
                        'N': [15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
                        'R': [0, 1, 0, 1, 0, 1],
                        'PHI': [10.0,
                                16.666666666666668,
                                16.666666666666668,
                                43.333333333333336,
                                43.333333333333336,
                                50.00000000000001]}})
    return femagtools.windings.Winding(wpar)


@ pytest.fixture()
def wdg1():
    wpar = dict(Q=12, p=5, m=3,
                windings={
                    1: {'N': [100.0, 100.0, 100.0, 100.0],
                        'layer': [1, 2, 1, 2], 'slots': [1, 1, -2, 6]},
                    2: {'N': [100.0, 100.0, 100.0, 100.0],
                        'layer': [2, 1, 2, 1], 'slots': [-4, 5, 5, -6]},
                    3: {'N': [100.0, 100.0, 100.0, 100.0],
                        'layer': [2, 1, 2, 1], 'slots': [2, -3, -3, 4]}})
    return femagtools.windings.Winding(wpar)


def test_slots(wdg):
    assert wdg.slots(1).tolist() == [
        [1,  1,  2,  5,  6,  6],
        [10, 10, 11, 14, 15, 15],
        [19, 19, 20, 23, 24, 24],
        [28, 28, 29, 32, 33, 33],
        [37, 37, 38, 41, 42, 42],
        [46, 46, 47, 50, 51, 51]]
    assert wdg.yd == 5
    assert wdg.zoneplan() == (
        [[1, 2, -6], [4, 5, -9], [-3, 7, 8]],
        [[1, -5, -6], [4, -8, -9], [-2, -3, 7]])


def test_axis(wdg):
    assert round(wdg.axis(), 3) == 0.349


def test_winding_factor(wdg):
    assert round(wdg.kw(), 4) == 0.9452


def test_winding_creation_1():
    wdg = femagtools.windings.Winding(dict(Q=12, p=2, m=3, l=1))
    assert wdg.slots(1).tolist() == [
        [1,  4], [7, 10]]
    assert wdg.zoneplan() == ([[1, -4], [3, -6], [-2, 5]], [])


def test_winding_creation_2():
    wdg = femagtools.windings.Winding(dict(Q=48, p=4, m=3, l=1))
    assert wdg.slots(1).tolist() == [
        [1,  2,  7,  8],
        [13, 14, 19, 20],
        [25, 26, 31, 32],
        [37, 38, 43, 44]]
    assert wdg.yd == 6


def test_custom_winding(wdg1):
    assert wdg1.slots(1).tolist() == [
        [1,  1,  2,  6,  7,  7,  8, 12]
    ]
    assert wdg1.zoneplan() == [
        [[1, 6, -7, -12], [-4, 5, 10, -11], [2, -3, -8, 9]],
        [[1, -2, -7, 8], [5, -6, -11, 12], [-3, 4, 9, -10]]]

    assert wdg1.yd == 1


def test_custom_winding_write(wdg1, tmp_path):
    name = 'winding'
    wdg1.write(name, tmp_path)
    p = tmp_path / (name + '.WID')
    assert 26 == len(p.read_text().split('\n'))


def test_coils_per_phase():
    wdg = femagtools.windings.Winding(
        {'Q': 48, 'p': 1, 'm': 3, 'l': 2, 'yd': 20})
    assert wdg.coils_per_phase() == 16


def test_turns_per_phase():
    wdg = femagtools.windings.Winding(
        {'Q': 48, 'p': 1, 'm': 3, 'l': 2, 'yd': 20})
    assert wdg.turns_per_phase(n=2, g=2) == 16


def test_inductance():
    wdg = femagtools.windings.Winding(
        {'Q': 36, 'p': 4, 'm': 3, 'l': 2, 'yd': 4})
    nc = 23
    g = 1
    lfe = 42e-3
    da1 = 110e-3
    ag = 1e-3
    assert round(wdg.inductance(nc, g, da1, lfe, ag), 4) == 0.0236


def test_unbalanced_winding():
    with pytest.raises(ValueError):
        wdg = femagtools.windings.Winding(
            {'Q': 48, 'p': 6, 'm': 3})
