from femagtools import amela
import json 
import os
import sys

def read_output(): 
    with open('tests/data/amela.out', 'r') as f: 
        data = f.readline()
    return data.split('\n')[0]

def read_json(): 
    with open('tests/data/pm_data/pm_data_se38.json', 'r') as f: 
        data = json.load(f)
    return data 

def test_amela():
    if sys.platform == 'linux':
        os.chmod('tests/data/AMELA', 0o0777)
    al = amela.Amela(workdir='tests/data', 
                    magnet_data=dict(name='pm_data'))
    loss = al()
    r = read_json()
    assert r['name'] == 'pm_data_se38'
    assert len(r['phi']) == 73  
    assert read_output() == 'succeed'
    assert round(loss['pm_data_se38']['total_loss'], 3) == 0.881                                             

    