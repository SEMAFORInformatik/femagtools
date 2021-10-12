"""
  Custom Winding definition demo
"""
machine = dict(
    name="SPM",
    lfe=0.1,
    poles=10,
    outer_diam=0.13,
    bore_diam=0.070,
    inner_diam=0.035,
    airgap=0.001,

    stator=dict(
        num_slots=12,
        statorRotor3=dict(
            slot_height=0.02,
            slot_h1=0.002,
            slot_h2=0.004,
            slot_r1=0.003,
            slot_r2=0.002,
            wedge_width1=0,
            wedge_width2=0,
            middle_line=1,  # vertical
            tooth_width=0.009,
            slot_top_sh=0.0,
            slot_width=0.003)
    ),

    magnet=dict(
        magnetSector=dict(
            magn_num=1,
            magn_height=0.004,
            magn_width_pct=0.8,
            condshaft_r=0.035,
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
        wdgdef={
            1: {'N': [100.0, 100.0, 100.0, 100.0],
                'layer': [1, 2, 1, 2], 'slots': [1, 1, -2, 6]},
            2: {'N': [100.0, 100.0, 100.0, 100.0],
                'layer': [2, 1, 2, 1], 'slots': [-4, 5, 5, -6]},
            3: {'N': [100.0, 100.0, 100.0, 100.0],
                'layer': [2, 1, 2, 1], 'slots': [2, -3, -3, 4]}})
    # num_phases=3,
    # num_wires=15,
    # coil_span=1,
    # num_layers=2)
)

if __name__ == '__main__':
    import femagtools
    import pathlib
    import logging
    import json

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    modelname = pathlib.Path(__file__).name.split('.')[0]
    logger = logging.getLogger(modelname)
    workdir = pathlib.Path().home() / 'femag'

    femag = femagtools.Femag(workdir)
    r = femag(machine)
    if r['status'] == 'ok':
        logger.info(f'{modelname} created')
    else:
        logger.error(f'{modelname} failed {r["message"]}')
