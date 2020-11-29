import femagtools


def create_fsl():
    machine = dict(
        name="statorfsl",
        desc="Stator FSL Demo",
        poles=10,
        outer_diam=0.1,
        bore_diam=0.055,
        inner_diam=0.033,
        airgap=0.001,
        lfe=0.15,
        stator=dict(
            num_slots=12,
            mcvkey_yoke='M270-35A',
            statorfsl=dict(
                yoke_height=0.008,
                slot_h1=0.0015,
                slot_h2=0.002,
                slot_width=0.003,
                tooth_width=0.007)
        ),
        magnet=dict(
            mcvkey_yoke="M270-35A",
            ring=dict(
                magn_height=0.005)
        ),
        windings=dict(
            num_phases=3,
            num_layers=2,
            num_wires=100,
            coil_span=1,
            cufilfact=0.4,
            culength=1.4)
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
    try:
        os.mkdir(workdir)
    except FileExistsError:
        pass

    with open(os.path.join(workdir, modelname+'.fsl'), 'w') as f:
        f.write('\n'.join(create_fsl()))

    logger.info("FSL %s created",
                os.path.join(workdir, modelname+'.fsl'))
