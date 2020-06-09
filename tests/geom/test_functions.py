import femagtools.dxfsl.geom as g
import numpy as np
import logging
import pytest


def test_alpha_points():
    p1 = (3.0, 1.0)
    p2 = (4.0, 5.0)
    p3 = (1.0, 3.0)

    alfa1 = g.alpha_points(p1, p2, p3)
    alfa2 = g.alpha_points(p2, p3, p1)
    alfa3 = g.alpha_points(p3, p1, p2)

    alfa = alfa1 + alfa2 + alfa3
    assert(np.isclose(alfa, 2*np.pi))


def test_points_are_close():
    rtol = 1e-05
    atol = 1e-08
    p1 = (3.14, 2.7666)
    p2 = (3.141592, 2.77)

    A = g.points_are_close(p1, p2, rtol, atol)
    B = g.points_are_close(p1, p2, rtol, 0.001)
    assert(A == B)

    A = g.points_are_close(p1, p2, 1e-01, atol)
    B = g.points_are_close(p1, p2, 1e-02, atol)
    assert(A == B)


def test_same_angle():
    assert(g.is_same_angle(3.14159266, np.pi))
    assert(g.is_same_angle(3.14159266, -np.pi))


def test_misc():
    assert(g.gcd(4, 5) == 1)
    assert(g.gcd(24, 6) == 6)
    assert(g.gcd(47, 23) == 1)
    assert(g.gcd(18, 3) == 3)
    assert(g.gcd(39, 13) == 13)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    test_alpha_points()
    test_points_are_close()
    test_same_angle()
    test_misc()
