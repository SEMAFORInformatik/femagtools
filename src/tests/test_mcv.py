import femagtools.mcv

def test_write_poc(tmpdir):
    filename = f"{tmpdir}/mcvData"

    B = [3e-02, 7e-02, 2e-01, 2.31]
    H = [1e+02, 2e+02, 5e+02, 3e+05]
    mcvData = dict(curve=[dict(bi=B,hi=H)])
    mcv = femagtools.mcv.MagnetizingCurve(mcvData)
    mcv.writefile(filename)
