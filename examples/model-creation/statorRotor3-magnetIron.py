import femagtools


def create_fsl():
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
            statorRotor3=dict(
                slot_height=0.02,
                slot_h1=0.002,
                slot_h2=0.004,
                slot_r1=0.0,
                slot_r2=0.0,
                wedge_width1=0,
                wedge_width2=0,
                middle_line=0,
                tooth_width=0.009,
                slot_top_sh=0.0,
                slot_width=0.003)
        ),

        magnet=dict(
            magnetIron=dict(
                magn_width=39e-3,
                magn_height=4e-3,
                gap_ma_iron=0.0,
                air_triangle=5e-3,
                iron_height=0.8e-3,
                magn_rem=1.2,
                condshaft_r=5e-3,
                bridge_height=0,
                magn_ori=1,
                bridge_width=0.0,
                iron_shape=33.5e-3)
        ),

        windings=dict(
            num_phases=3,
            num_wires=100,
            coil_span=3.0,
            num_layers=1)
    )

    return femagtools.create_fsl(machine)
    

if __name__ == '__main__':
    import os
    import logging
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    modelname = os.path.split(__file__)[-1].split('.')[0]
    logger = logging.getLogger(modelname)
    workdir = os.path.join(os.path.expanduser('~'), 'femag')

    with open(os.path.join(workdir, modelname+'.fsl'), 'w') as f:
        f.write('\n'.join(create_fsl()))

    logger.info("FSL %s created",
                os.path.join(workdir, modelname+'.fsl'))
