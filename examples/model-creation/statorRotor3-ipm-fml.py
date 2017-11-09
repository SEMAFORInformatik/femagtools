import os
import femagtools

machine = dict(
    name="IPM 270 L8",
    desc="IPM FML Demo",
    poles=8,
    outer_diam=0.26924,
    bore_diam=0.16192,
    inner_diam=0.11064,
    airgap=0.00075,
    lfe=0.08356,
    stator=dict(
        num_slots=48,
        mcvkey_yoke='M270-35A',
        zeroangle=67.5,
        statorRotor3=dict(
            slot_height=0.0335,
            slot_h1=0.001,
            slot_h2=0.0,
            slot_width=0.00193,
            slot_r1=0.0001,
            slot_r2=0.00282,
            wedge_width1=0.00295,
            wedge_width2=0.0,
            middle_line=0.0,
            tooth_width=0.0,
            slot_top_sh=0.0)
    ),
    magnet=dict(
        mcvkey_yoke="M270-35A",
        magnetFsl=dict(
            magn_height=0.004,
            magn_width=0.042,
            gap_ma_iron=0.009,
            gap_air_iron=0.0015,
            gamma=49,
            content_template="ipm-fml.fsl")
    ),
    windings=dict(
        num_phases=3,
        num_layers=1,
        num_wires=9,
        coil_span=6.0,
        cufilfact=0.4,
        culength=1.4,
        slot_indul=0.5e-3
    )
)

workdir = os.path.join(os.path.expanduser('~'), 'femag')
with open(os.path.join(workdir, 'stator3-ipm.fsl'), 'w') as f:
    f.write('\n'.join(femagtools.create_fsl(machine)))
