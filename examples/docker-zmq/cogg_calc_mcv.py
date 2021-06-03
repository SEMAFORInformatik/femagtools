import femagtools
import os
import logging

FORMAT = "%(asctime)s %(name)s: %(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger("main")


def show_progress(message):
    """handle messages from femag publisher"""
    topic, content = message  # "femag_log" + text
    if topic == 'femag_log' or topic == 'progress':
        logger.info("%s: %s",
                    topic, content.strip())


mcvData = [
    dict(curve=[dict(
        bi=[0.0, 0.09, 0.179, 0.267, 0.358,
            0.45, 0.543, 0.6334, 0.727,
            0.819, 0.9142, 1.0142, 1.102,
            1.196, 1.314, 1.3845, 1.433,
            1.576, 1.677, 1.745, 1.787,
            1.81, 1.825, 1.836],

        hi=[0.0, 22.16, 31.07, 37.25, 43.174,
            49.54, 56.96, 66.11, 78.291,
            95, 120.64, 164.6, 259.36,
            565.86, 1650.26, 3631.12, 5000, 10000,
            15000, 20000, 25000, 30000, 35000, 40000])],
         name='m270-35a',
         desc=u"Demo Steel",
         ch=0.0,
         cw_freq=2.0,
         cw=1.68),

    dict(curve=[dict(
        bi=[-0.0178, 0.08683, 0.19248, 0.29856, 0.40417,
            0.51, 0.6157, 0.72186, 0.8281, 0.93357, 1.0396,
            1.14526, 1.240],
        
        hi=[-878813.5,
            -848291.125,
            -804114.125,
            -747085.625,
            -681221.6875,
            -610538.375,
            -536642.1875,
            -461139.59375,
            -384030.53125,
            -306118.21875,
            -228205.9375,
            -149490.421875,
            -70774.90625])],
         name='BM38H',
         remz=1.24,
         ctype=2)]

mcv = femagtools.mcv.MagnetizingCurve(mcvData)

magnetmat = [
    dict(name='BM38H',
         mcvkey='BM38H',
         orient='mpolaniso')]

machine = dict(
    name="PM 225 8",
    lfe=0.1,
    poles=8,
    outer_diam=0.225,
    bore_diam=0.1615,
    inner_diam=0.120,
    airgap=0.0015,
     
    stator=dict(
        num_slots=48,
        mcvkey_yoke="m270-35a",
        rlength=1.0,
        statorRotor3=dict(
            slot_height=0.0197,
            slot_h1=0.001,
            slot_h2=0.001,
            slot_width=0.003,
            slot_r1=0.0,
            slot_r2=0.0,
            wedge_width1=0.0,
            wedge_width2=0.0,
            tooth_width=0.005,
            middle_line=1,
            slot_top_sh=1)
    ),
    
    magnet=dict(
        mcvkey_shaft="dummy",
        material='BM38H',
        mcvkey_yoke="m270-35a",
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
            magn_len=1.0)
    ),

    windings=dict(
        num_phases=3,
        num_wires=4,
        coil_span=5,
        slot_indul=0,
        num_layers=2)
    )


workdir = os.path.join(
    os.path.expanduser('~'), 'femag')
try:
    os.makedirs(workdir)
except OSError:
    pass

femag = femagtools.ZmqFemag(5555, workdir=workdir,
                            magnets=magnetmat,
                            magnetizingCurves=mcv)

simulation = dict(
    calculationMode="cogg_calc",
    magn_temp=60.0,
    period_frac=6,
    speed=50.0)

r = femag(machine, simulation, show_progress)

print("Order    T/Nm      %")
tq = r.torque_fft[-1]
for l in zip(tq['order'], tq['torque'], tq['torque_perc']):
    if(l[2] > 1):
        print('{0:<5} {1:9.2f} {2:6.1f}'.format(*l))
# remove temp files
femag.cleanup()
