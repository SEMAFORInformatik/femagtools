import sys
import os
import pytest
import numpy as np
import femagtools.mcv
import pathlib

@pytest.fixture
def mcv_data_dir():
    return pathlib.Path(__file__).with_name('data')

def test_write_mcv_magcrv(tmpdir):
    filename = f"{tmpdir}/mcvData"

    B = [3e-02, 7e-02, 2e-01, 2.31]
    H = [1e+02, 2e+02, 5e+02, 3e+05]
    mcvData = dict(curve=[dict(bi=B,hi=H)])
    mcv = femagtools.mcv.MagnetizingCurve(mcvData)
    mcv.writefile(filename)
    mcvr = femagtools.mcv.Reader()
    ext = '.MC' if sys.platform == 'win32' else '.MCV'
    mcvr.readMcv(filename+ext)
    assert mcvr['version_mc_curve'] == 0
    mcvr['ctype'] == femagtools.mcv.MAGCRV
    assert mcvr['curve'][0]['bi'] == pytest.approx(B, 1e-3)
    assert mcvr['curve'][0]['hi'] == pytest.approx(H, 1e-3)
    assert len(mcvr['curve']) == 1

def test_write_mcv_orient_crv(tmpdir):
    filename = f"{tmpdir}/mcvData"

    B = [3e-02, 7e-02, 2e-01, 2.31]
    H = [1e+02, 2e+02, 5e+02, 3e+05]
    mcvData = dict(curve=[dict(bi=B,hi=H), dict(bi=B,hi=H)])
    mcv = femagtools.mcv.MagnetizingCurve(mcvData)
    mcv.writefile(filename)
    ext = '.MC' if sys.platform == 'win32' else '.MCV'
    mcvr = femagtools.mcv.Reader()
    mcvr.readMcv(filename+ext)
    assert mcvr['version_mc_curve'] == 1
    assert mcvr['ctype'] == femagtools.mcv.MAGCRV
    assert mcvr['curve'][0]['bi'] == pytest.approx(B, 1e-3)
    assert mcvr['curve'][0]['hi'] == pytest.approx(H, 1e-3)
    assert mcvr['curve'][1]['bi'] == pytest.approx(B, 1e-3)
    assert mcvr['curve'][1]['hi'] == pytest.approx(H, 1e-3)
    assert len(mcvr['curve']) == 2

def test_read_mcv_magcrv(mcv_data_dir):
    filename = "TKS_NO_20.MCV"
    reader = femagtools.mcv.Reader()
    reader.readMcv(mcv_data_dir / filename)
    r = reader.get_results()
    assert r['desc'] == u'PowerCore\xae NO 20 ;ThyssenKrupp Steel Eur'
    assert r['fillfac'] == pytest.approx(0.92, 1e-6)
    assert len(r['curve'][0]['bi']) == 24
    assert r['curve'][0]['bi'][0] == 0.0
    assert r['curve'][0]['bi'][-1] == pytest.approx(1.836, 1e-3)

def test_read_write_mcv_with_losses(mcv_data_dir):
    filename = "TKM270-50A-LOSS.MCV"
    reader = femagtools.mcv.Reader()
    reader.readMcv(mcv_data_dir / filename)
    r = reader.get_results()
    assert r['name'] == u'TKM270-50A-LOSS'
    assert len(r['curve'][0]['bi']) == 35
    assert r['curve'][0]['bi'][-1] == pytest.approx(2.638, 1e-3)
    assert set(r['losses'].keys()) == {'B', 'f', 'pfe', 'Bo', 'fo',
                                       'cw', 'cw_freq', 'b_coeff', 'ch', 'ch_freq'}
    assert len(r['losses']['pfe']) == 20
    assert len(r['losses']['pfe'][8]) == 19
    assert r['losses']['pfe'][8][18] == pytest.approx(3097.6, 0.1)

    # test mcv writer
    filename_out = "TKS_LOSS.MCV"
    writer = femagtools.mcv.Writer(r)
    writeMcvFile = mcv_data_dir / filename_out
    writer.writeMcv(writeMcvFile)
    # TEST
    reader2 = femagtools.mcv.Reader()
    reader2.readMcv(writeMcvFile)

    for attr in ['version_mc_curve', 'mc1_curves', 'mc1_title']:
        assert getattr(reader, attr) == getattr(reader2, attr)
    for attr in ['mc1_remz', 'mc1_bsat', 'mc1_bref', 'mc1_fillfac',
                 'fo', 'Bo', 'ch', 'ch_freq', 'cw',
                 'cw_freq', 'b_coeff', 'rho', 'fe_sat_mag']:
        assert getattr(reader, attr) == getattr(reader2, attr)
    for attr in ['hi', 'bi']:
        assert reader.curve[0][attr] == reader2.curve[0][attr]
    np.testing.assert_almost_equal(reader.losses['f'],
                                   reader2.losses['f'], 5)
    np.testing.assert_almost_equal(reader.losses['B'],
                                   reader2.losses['B'], 5)

    writer.bertotti = {'ch': 0.02196930496538005,
                      'cw': 9.739048407022554e-05,
                       'ce': 0.00047834510925596726,
                       'ch_freq': 2, 'b_coeff':0}
    writer.writeMcv(mcv_data_dir / "TKS_LOSS2.MCV", feloss='bertotti')


def test_extrapol():
    mcv = {'name':"M222",
           'curve':[{
            'bi':[0.0, 0.09, 0.179, 0.267, 0.358,
                0.45, 0.543, 0.6334, 0.727,
                0.819, 0.9142, 1.0142, 1.102,
                1.196, 1.314, 1.3845, 1.433,
                1.576, 1.677, 1.745, 1.787,
                1.81, 1.825, 1.836],

            'hi':[0.0, 22.16, 31.07, 37.25, 43.174,
                49.54, 56.96, 66.11, 78.291,
                95, 120.64, 164.6, 259.36,
                565.86, 1650.26, 3631.12, 5000, 10000,
                15000, 20000, 25000, 30000, 35000, 40000]}]}
    c = femagtools.mcv.extrapol(mcv['curve'][0])
    assert len(c['bi']) > len(mcv['curve'][0]['bi'])
    assert len(c['hi']) > len(mcv['curve'][0]['hi'])
    assert np.all(np.diff(c['bi']) > 0)
    assert np.all(np.diff(c['hi']) > 0)
