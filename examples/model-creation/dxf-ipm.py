import femagtools

def create_fsl():
    machine = dict(
        name="IPM4",
        desc="DXF Demo",
        lfe=0.08356,
        
        dxffile=dict(
            name='ipm4.dxf'
        ),

        stator=dict(
            mcvkey_yoke='M270-35A'
        ),
        magnet=dict(
            mcvkey_yoke="M270-35A"
        ),
    
        windings=dict(
            num_phases=3,
            num_layers=2,
            num_wires=10,
            coil_span=8
        )
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
