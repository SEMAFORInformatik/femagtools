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
            num_slots_gen=3,
            mcvkey_yoke="dummy",
            rlength=1.0,
            stator1=dict(
                slot_rf1=0.057,
                tip_rh1=0.037,
                tip_rh2=0.037,
                tooth_width=0.009,
                slot_width=0.003)
        ),

        magnet=dict(
            mcvkey_shaft="dummy",
            mcvkey_yoke="dummy",
            magnetSector=dict(
                magn_num=1,
                magn_width_pct=0.8,
                magn_height=0.004,
                magn_shape=0.0,
                bridge_height=0.0,
                magn_type=1,
                condshaft_r=0.02,
                magn_ori=2,
                magn_rfe=0.0,
                bridge_width=0.0,
                magn_len=1.0)
        ),

        windings=dict(
            num_phases=3,
            num_wires=100,
            coil_span=3.0,
            num_layers=1)
    )

    operatingConditions = dict(
        calculationMode="pm_sym_fast",
        current=12.0,
        angl_i_up=0.0,
        speed=15.0,
        wind_temp=60.0,
        magn_temp=60.0)

    return femagtools.create_fsl(machine, operatingConditions)


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

