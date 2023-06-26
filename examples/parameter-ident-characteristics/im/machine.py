"""
Parameter Identification of the equivalent electrical circuit
for an Induction Machine based on Femag AC and femagtools

Expected values:
 pnom = 5e3  -- nominal power W
 f1nom = 50.0  -- Frequency Hz
 speed = 1450 -- speed 1/min  (s=3.33%)
 u1nom = 400  -- nominal line voltage V
 i1nom = 10.4 -- nominal line current A
  i10  = 6.7 -- no load current A
  conn = star
  Tnom = 32.9 -- Nm

Source: Pavel Ponomarev (May 2017)
  https://www.researchgate.net/publication/317012206_SEMTEC_Report_Elmer_FEM_-_Induction_Machine_Tutorial

Ronald Tanner
"""

import femagtools.mcv
#from femagtools.docker import Engine


magnetizingCurves = femagtools.mcv.MagnetizingCurve([
    {"name": "M800-65A",
     "desc": "for demo purposes only",
     "ctype": 1, "fillfac": 1.0, "fo": 50.0, "Bo": 1.5,
     "ch": 0.0, "ch_freq": 1,  "cw": 3, "cw_freq": 1.45, "b_coeff": 1.0,
     "rho": 7.65,
     "curve": [
         {"hi": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0,
                 80.0, 90.0, 100.0, 200.0, 300.0, 400.0, 500.0,
                 600.0, 700.0, 800.0, 900.0, 1000.0, 2000.0,
                 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0,
                 9000.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0,
                 60000.0, 70000.0, 80000.0, 90000.0, 100000.0, 200000.0,
                 500000.0],
          "bi": [0.0, 0.056587, 0.112851, 0.168478, 0.223178, 0.2766881, 0.328783,
                 0.379275, 0.428017, 0.474904, 0.519865, 0.866555, 1.072363,
                 1.2, 1.285121, 1.345424, 1.39036, 1.4252, 1.453078, 1.475977,
                 1.591754, 1.644, 1.68, 1.71, 1.735142, 1.758716, 1.78075, 1.80158,
                 1.82141, 1.9778, 2.0735, 2.129, 2.162473, 2.184886, 2.2021,
                 2.2169, 2.2306, 2.243757, 2.37071, 2.7482841],
          "angle": 0.0}]
     }])

condMat = [
    {"name": "Al", "spmaweight": 2.7, "elconduct": 33e6, "tempcoef": 3.9e-3},
    {"name": "Cu", "spmaweight": 8.96, "elconduct": 56e6, "tempcoef": 3.9e-3}
]

machine = dict(
    name="VTT-PP",
    f1nom=50, u1nom=230,
    lfe=0.16,
    poles=4,
    outer_diam=0.220,
    bore_diam=0.125,
    shaft_diam=0.056,
    inner_diam=0.005,
    # inner_diam=0.056,
    airgap=0.5e-3,

    stator=dict(
        num_slots=48,
        mcvkey_yoke='M800-65A',
        fillfac=0.97,
        nodedist=2.0,
        statorRotor3=dict(
            slot_height=0.0283,
            slot_h1=7e-4,
            slot_h2=9e-4,
            slot_r1=0.0,
            slot_r2=3.4e-3,
            wedge_width1=4e-3,
            wedge_width2=0,
            tooth_width=0.0,
            slot_top_sh=0.0,
            slot_width=2.8e-3
        )
    ),

    rotor=dict(
        num_slots=40,
        material='Al',
        mcvkey_yoke='M800-65A',
        fillfac=0.97,
        conductivity=33.4e6,  # Al in S/m
        statorRotor3=dict(
            slot_height=0.0156,
            slot_h1=0.5e-3,
            slot_h2=2.2e-3,
            slot_r1=2.2e-3,
            slot_r2=1e-3,
            wedge_width1=1.0e-3,
            wedge_width2=0,
            tooth_width=0.0,
            slot_top_sh=0.0,
            slot_width=1e-3
        )

    ),

    windings=dict(
        material='Cu',
        num_par_wdgs=2,
        dia_wire=1.5e-3,
        num_phases=3,
        num_layers=1,
        num_wires=32,
        coil_span=12,
        fillfac=0.42
    )
)
