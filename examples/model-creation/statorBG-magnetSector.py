import femagtools

def create_fsl():
    machine = dict(
        name="PM 130 L4",
        lfe=0.1,
        poles=4,
        outer_diam=0.04,
        bore_diam=0.022,
        inner_diam=0.005,
        airgap=0.001,

        stator=dict(
            num_slots=12,
            statorBG=dict(
                yoke_diam_ins=0.0344,
                slot_h1=0.0005,
                slot_h3=0.004,
                slot_width=0.00395,
                slot_r1=0.0008,
                slot_r2=0.0007,
                tooth_width=0.0032,
                middle_line=1,
                tip_rad=-0.03,
                slottooth=0)
        ),

        magnet=dict(
            magnetSector=dict(
                magn_num=1,
                magn_width_pct=0.8,
                magn_height=0.002,
                magn_shape=0.0,
                bridge_height=0.0,
                magn_type=1,
                condshaft_r=0.005,
                magn_ori=2,
                magn_rfe=0.0,
                bridge_width=0.0,
                magn_len=1.0)
        ),

        windings=dict(
            num_phases=3,
            num_wires=10,
            coil_span=1,
            num_layers=2)
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
