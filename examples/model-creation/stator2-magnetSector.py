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
            num_slots_gen=12,
            stator2=dict(
                slot_t1=0.001,
                slot_t2=0.0005,
                slot_t3=0.0005,
                slot_depth=0.020,
                slot_width=0.009,
                corner_width=0.01)
        ),

        magnet=dict(
            magnetSector=dict(
                magn_num=1,
                magn_width_pct=0.8,
                magn_height=0.004,
                magn_shape=0.0,
                bridge_height=0.0,
                magn_type=1,
                condshaft_r=0.02,
                magn_ori=8,
                magn_rfe=0.0,
                bridge_width=0.0,
                magn_len=1.0)
        ),

        winding=dict(
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
