import femagtools.heat_source_network as hsn
import pathlib
import pytest


@pytest.fixture
def data_dir():
    return pathlib.Path(__file__).with_name('data')


def test_heat_source_network(data_dir):
    model = hsn.read(str(data_dir / 'temp_model.hsn'))
    assert [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13] == model.nodes
    names = model.get_node_names()
    assert ['StZa', 'outs', 'StJo', 'Slot',
            'Shaf', 'Iron', 'PMag', 'PMag',
            'W1  ', 'W2  ', 'W3  '] == names
    T = model.solve()
    # [130.22552926,  56.75737966, 116.12788806, 150.18881253,
    #   175.77455152, 175.77455152, 175.77526274, 175.77455152,
    #   126.04691063, 122.49609205, 125.27401647])
