import femagtools

machine = dict(
    name="DASSW",
    description="Dual Airgap, Single Stator Winding",
    lfe=0.023,
    poles=4,
    outer_diam=0.072,
    bore_diam=[0.025, 0.0558],
    inner_diam=0.015,
    airgap=[0.001,0.001],
    external_rotor=True,

    stator=dict(
        num_slots=4,
        stator=dict(
            beta_s=0.3,
            nodepitch1=1.0,
            nodepitch2=0.5)
    ),
    magnet=dict(
        remancenc=0.415,
        relperm=1.09,
        magnet=dict(
            hm=5e-3,
            alpha_m=0.98)
    ),
    windings=dict(
        num_phases=1,
        num_layers=2,
        num_wires=500,
        coil_span=1)
)

with open(machine['name'] +'.fsl', 'w') as fp:
    fp.write('\n'.join(femagtools.create_fsl(machine)))
