import femagtools


def create_fsl():
    machine = dict(
        name="example2",
        lfe=0.1,
        poles=4,
        outer_diam=0.143,
        bore_diam=0.07,
        inner_diam=0.015,
        airgap=0.001,

        stator=dict(
            num_slots=12,
            mcvkey_yoke="dummy",
            rlength=1.0,
            stator1=dict(
                slot_rf1=0.057,
                tip_rh1=0.037,
                tip_rh2=0.037,
                tooth_width=0.011,
                slot_width=0.003)
        ),

        magnet=dict(
            mcvkey_yoke="dummy",
            spokefml=dict(
                magn_height=0.008,
                shaft_diam=0.01,
                slot_width=0.004,
                magn_width=0.024
            )
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
