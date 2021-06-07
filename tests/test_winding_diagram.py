
from femagtools import winding_diagram


def test_winding_diagram():

    data = winding_diagram._winding_data(12, 2, 3)
    assert data == [1, -2, 3, -1, 2, -3, 1, -2, 3, -1, 2, -3]

    data = winding_diagram._winding_data(36, 2, 3)
    assert data == [1, 1, 1, -2, -2, -2, 3, 3, 3, -1, -1, -1, 2, 2, 2, -3, -3, -3, 1, 1, 1, -2, -2, -2, 3, 3, 3, -1, -1, -1, 2, 2, 2, -3, -3, -3]

    data = winding_diagram._winding_data(36, 3, 3)
    assert data == [1, 1, -2, -2, 3, 3, -1, -1, 2, 2, -3, -3, 1, 1, -2, -2, 3, 3, -1, -1, 2, 2, -3, -3, 1, 1, -2, -2, 3, 3, -1, -1, 2, 2, -3, -3]

