"""Electrical Excited Synchronous Machine Example

  Author: Ronald Tanner, Semafor Informatik & Energie AG, Switzerland
"""
# DISS. ETH NO. 22393
# Illiano, Enzo Michele
# DESIGN OF A HIGHLY EFFICIENT BRUSHLESS CURRENT
# EXCITED SYNCHRONOUS MOTOR FOR AUTOMOTIVE PURPOSES
#                 A     B    C     D     E      F
# Iq (A)         231    45   231   380   185     65
# Id (A)         -28   -50   -28   -90  -130   -130
# Iex (A)        3.2   1.1   3.2   5.4   2.8    3.2
# I1            164.5  47.6 164.5  276.1  159.9 102.8
# Beta °        -7     -48   -7    -13.2 -35.1  -63.4
# Speed  (rpm)  2000  4000  4000  4000  6000  10000
# Torque (Nm)    120    20   120   220   100     60
# VOLTAGE, (V)  93.1  77.0  183.5 198.1 213.5 189.9

import pathlib
import json

laminations = [
    json.load(p.open()) for p in pathlib.Path(
        './lamination').glob('*.json')]


condMat = [
    {"name": "Al", "spmaweight": 2.7, "elconduct": 33e6, "tempcoef": 3.9e-3},
    {"name": "Cu", "spmaweight": 8.96, "elconduct": 56e6, "tempcoef": 3.9e-3}
]

machine = dict(
    name="EESM",
    lfe=0.123,
    poles=6,
    outer_diam=0.240,
    bore_diam=0.1658,
    inner_diam=0.085,
    airgap=5e-4,
    # num_agnodes=960 --> m.delta_angle_ndchn = 360/num_agnodes --   angle nodechain airgap in [degr]

    stator=dict(
        mcvkey_yoke='M270-35A',
        nodedist=2,
        num_slots=54,
        statorRotor3=dict(
            slot_height=0.02065,
            slot_h1=0.001,
            slot_h2=0.001,
            slot_r1=0.0,
            slot_r2=0.0,
            wedge_width1=0,
            wedge_width2=0,
            middle_line=0,
            tooth_width=0.0057,
            slot_top_sh=0.0,
            slot_width=0.0024)
    ),

    rotor=dict(
        mcvkey_yoke='M270-35A',
        num_wires=400,
        resistance=40,
        nodedist=2,
        ifnom=6.0,
        rot_hsm=dict(
            gap_pol_shaft=0.0,
            core_height=0.0165,
            pole_height=0.0115,
            pole_rad=0.0725,
            core_width2=0.0322,
            core_width1=0.0322,
            pole_width_r=0.076,
            pole_width=0.075,
            slot_width=0,
            slot_height=0,
            damper_diam=0,
            damper_div=0
        )
    ),

    windings=dict(
        material='Cu',
        num_par_wdgs=3,
        num_phases=3,
        num_wires=7,
        coil_span=9,
        num_layers=1)
)
