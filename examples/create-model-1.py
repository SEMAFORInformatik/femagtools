import femagtools

machine = dict(
    name="PM 886 32",
    lfe=0.224,
    poles=32,
    outer_diam=0.886,
    bore_diam=0.76,
    inner_diam=0.4956,
    airgap=0.007,
    external_rotor=1,
     
    stator=dict(
        num_slots=120,
        mcvkey_yoke="dummy",
        rlength=1.0,
        stator4=dict(
            slot_height=0.035,
            slot_h1=0.002,
            slot_h2=0.0,
            slot_h3=0.004,
            slot_h4=0.0,
            slot_width=0.01,
            slot_r1=0.0,
            wedge_width1=0.01,
            wedge_width2=0.01,
            wedge_width3=0.01)
    ),
    
    magnet=dict(
        mcvkey_yoke="dummy",
        magnetSector=dict(
            magn_num=1,
            magn_height=0.014,
            magn_width_pct=0.85,
            condshaft_r=0.0,
            magn_rfe=0.0,
            magn_len=1.0,
            magn_shape=0.0,
            bridge_height=0.0,
            bridge_width=0.0,
            magn_ori=1,
            magn_type=1
        )
    ),

    windings=dict(
        num_phases=3,
        num_wires=5,
        coil_span=1,
        num_layers=2)
)

fsl = femagtools.create_fsl(machine)
print('\n'.join(fsl))
