machine = {
    'name': "PM 130 L4",
    'lfe': 0.05,
    'poles': 4,
    'outer_diam': 0.11,
    'bore_diam': 0.057,
    'inner_diam': 0.015,
    'airgap': 0.0015,
    'stator': {
        'num_slots': 6,
        'mcvkey_yoke': 'M400_50A',
        'rlength': 1,
        'simplestator': {
            'slot_height': 16e-3,
            'slot_width': 16e-3,
        },
    },
    'magnet': {
        'mcvkey_yoke': 'Stahl37',
        'rlength': 1,
        'simplemagnet': {
            'magnet_height': 5e-3,
            'magnet_width': 30e-3
        }
    },
    'winding': {
        'num_phases': 3,
        'num_wires': 10,
        'coil_span': 3,
        'num_layers': 1
    }
}
