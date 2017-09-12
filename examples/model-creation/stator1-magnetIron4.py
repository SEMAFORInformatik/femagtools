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
        magnetIron4=dict(
            magn_height=4e-3,
            magn_width=36e-3,
            gap_ma_iron=1e-3,
            air_space_h=3e-3,
            iron_bfe=2.5e-3,
            magn_di_ra=8e-3,
            corner_r=2e-3,
            air_sp_ori=0.0,
            condshaft_r=7e-3,
            magn_ori=1,
            magn_num=2,
            iron_shape=0)
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

