import femagtools.femag


def test_run_script(monkeypatch, tmpdir):
    def mock_run(*args, **kwargs):
        return
    monkeypatch.setattr(femagtools.femag.Femag, "run", mock_run)
    femag = femagtools.femag.Femag(tmpdir)
    r = femag(dict(), dict())
    assert r['status'] == 'ok'
  
