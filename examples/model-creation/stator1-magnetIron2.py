import femagtools

machine = dict(
    name="PM 130 L4",
    lfe=0.1,
    poles=4,
    outer_diam=0.13,
    bore_diam=0.07,
    inner_diam=0.015,
    airgap=0.001,
     
    stator=dict(
        num_slots=12,
        stator1=dict(
            slot_rf1=0.057,
            tip_rh1=0.037,
            tip_rh2=0.037,
            tooth_width=0.009,
            slot_width=0.003)
    ),
    
    magnet=dict(
        magnetIron2=dict(
            magn_width=39e-3,
            magn_height=4e-3,
            gap_ma_iron=1e-3,
            air_triangle=1e-3,
            iron_height=1e-3,
            gap_ma_right=0.0,
            gap_ma_left=0.0,
            magn_rem=1.2,
            condshaft_r=5e-3,
            magn_ori=1,
            iron_shape=33.5e-3)
    ),
    
    windings=dict(
        num_phases=3,
        num_wires=100,
        coil_span=3.0,
        num_layers=1)
)

fsl = femagtools.create_fsl(machine)
with open('femag.fsl', 'w') as f:
    f.write('\n'.join(fsl))

