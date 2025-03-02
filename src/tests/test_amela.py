import pytest
import json
from femagtools import amela
from pathlib import Path

def read_output():
    return (Path(__file__).parent/ "data/amela.out").read_text().split('\n')[0]

def read_json():
    return json.loads(
        (Path(__file__).parent/"data/pm_data/pm_data_se38.json").read_text())

def test_amela():
    al = amela.Amela(workdir='src/tests/data',
                     magnet_data=dict(name='pm_data'))
    loss = al()
    r = read_json()
    assert r['name'] == 'pm_data_se38'
    assert len(r['phi']) == 73
    assert read_output() == 'succeed'
    assert round(loss['pm_data_se38']['total_loss'], 3) == 0.881
