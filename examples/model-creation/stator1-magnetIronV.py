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
            stator1=dict(
                slot_rf1=0.057,
                tip_rh1=0.037,
                tip_rh2=0.037,
                tooth_width=0.009,
                slot_width=0.003)
        ),

        magnet=dict(
            magnetIronV=dict(
                magn_height=4e-3,
                magn_width=18e-3,
                magn_angle=130,
                iron_hs=1e-3,
                iron_height=2e-3,
                gap_ma_iron=1e-3,
                air_triangle=2e-3,
                condshaft_r=12e-3,
                magn_rem=1.2,
                magn_num=2,
                iron_shape=0)
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

