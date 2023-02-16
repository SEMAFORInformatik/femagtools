import femagtools.femag
import femagtools.mcv

mcvData = [
    dict(curve=[dict(
        bi=[0.0, 0.45, 0.819, 1.196, 1.576, 1.81],
        hi=[0.0, 49.54, 95, 565.86, 15000, 20000])],
         name='m270-35a',
         desc=u"Demo Steel",
         ch=4.0,
         cw_freq=2.0,
         cw=1.68)]
mcv = femagtools.mcv.MagnetizingCurve(mcvData)

magnetmat = [
    dict(
        name='M45',
        remanenc=1.1,
        relperm=1.04,
        spmaweight=7.4,
        temcoefbr=-0.0015,
        temcoefhc=-0.0013,
        magncond=625000.0
    )]

machine = dict(
    name="PM 130 L4",
    outer_diam=0.13,
    bore_diam=0.07,
    inner_diam=0.,
    lfe=0.1,
    poles=4,
    stator=dict(
        mcvkey_yoke='m270-35a',
        num_slots=12,
        num_slots_gen=3,
        nodedist=1.5,
        rlength=1.0),
    magnet=dict(
        material='M45',
        magnetSector=dict(
            magn_num=1,
            magn_width_pct=0.8,
            magn_height=0.005,
            magn_shape=0.0,
            bridge_height=0.0,
            magn_type=1,
            condshaft_r=0.02,
            magn_ori=2,
            magn_rfe=0.0,
            bridge_width=0.0,
            magn_len=1.0)),
    windings=dict(
        num_phases=3,
        num_layers=1,
        num_wires=4,
        coil_span=3))

simulation = dict(
    calculationMode="cogg_calc",
    magn_temp=60.0,
    speed=50.0)


def test_run_script(monkeypatch, tmpdir):
    def mock_run(*args, **kwargs):
        return
    monkeypatch.setattr(femagtools.femag.Femag, "run", mock_run)
    femag = femagtools.femag.Femag(str(tmpdir),
                                   magnetizingCurves=mcv, magnets=magnetmat)

    fsl = femag.create_fsl(machine, simulation)
    assert len(fsl) == 92
    assert femag.model.magnet['temp_prop']['magntemp'] == 60.0

    r = femag(machine, simulation)

    assert isinstance(r, femagtools.bch.Reader)
    assert tmpdir.join("femag.fsl").exists()
