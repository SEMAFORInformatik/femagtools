#!/usr/bin/env python
#
import pytest
import numpy as np
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
    assert wdg.yd == 4
    assert wdg.zoneplan() == (
        [[1, 2, -6], [4, 5, -9], [-3, 7, 8]],
        [[1, -5, -6], [4, -8, -9], [-2, -3, 7]])


def test_axis():
    for wpar, expected in [
            ({'Q': 168, 'p': 7, 'm': 3, 'l': 2, 'yd': 10}, 14.99973842075893),
            ({'Q': 90, 'p': 12, 'm': 3, 'l': 2, 'yd': 4},  11.99951171875),  # TODO
            ({'Q': 108, 'p': 9, 'm': 3, 'l': 2, 'yd': 5},  11.666259765625),
            ({'Q': 96, 'p': 4, 'm': 3, 'l': 2, 'yd': 10},  26.249542236328125),
            ({'Q': 144, 'p': 16, 'm': 3, 'l': 2, 'yd': 4},  7.49969482421875),
            ({'Q': 144, 'p': 12, 'm': 3, 'l': 2, 'yd': 5},  8.74969482421875),
            ({'Q': 54, 'p': 3, 'm': 3, 'l': 2, 'yd': 8},  36.66585286458334),
            ({'Q': 48, 'p': 4, 'm': 3, 'l': 1, 'yd': 6}, 29.99908447265625)]:
        wdg = femagtools.windings.Winding(wpar)
        assert round(wdg.axis()/np.pi*180, 2) == round(expected, 2)


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


def test_num_layers():
    w = femagtools.windings.Winding(
        {"Q": 48, "p": 4, "m": 3, "windings": {
            1: {"N": [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                "layer": [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
                "slots": [-5, 1, -6, 2, 11, -7, 12, -8]},
            2: {"N": [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                "layer": [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
                "slots": [3, 5, 4, 6, -9, -11, -10, -12]},
            3: {"N": [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
                "layer": [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2],
                "slots": [-1, -3, -2, -4, 7, 9, 8, 10]}}})

    assert w.l == 2
    assert w.yd == 4

def test_zoneplan():
    w = femagtools.windings.Winding(
        {'Q': 60, 'p': 32, 'm': 3, 'l': 1})
    assert [1, -2, 3, -4, 5, -6, -17, 18, -19, 20] == w.zoneplan()[0][0]
