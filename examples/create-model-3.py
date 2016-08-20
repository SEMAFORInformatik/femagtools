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
        num_slots_gen=3,
        mcvkey_yoke="dummy",
        rlength=1.0,
        stator1=dict(
            slot_rf1=0.057,
            tip_rh1=0.037,
            tip_rh2=0.037,
            tooth_width=0.009,
            slot_width=0.003)
    ),
    
    magnet=dict(
        mcvkey_yoke="dummy",
        magnetFsl=dict(
            magn_height=0.014,
            content_template="ring-magnet.fsl"
        )
    ),
    
    windings=dict(
        num_phases=3,
        num_wires=100,
        coil_span=3.0,
        num_layers=1)
)

operatingConditions = dict(
    calculationMode="pm_sym_fast",
    current=50.0,
    angl_i_up=0.0,
    speed=50.0,
    wind_temp=60.0,
    magn_temp=60.0)

fsl = femagtools.create_fsl(machine, operatingConditions)
print('\n'.join(fsl))
