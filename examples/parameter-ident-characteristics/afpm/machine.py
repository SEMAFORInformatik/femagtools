""" Axial Flux PM Machine (AFPM) example

 Suppported types
 "S1R1"      -- 1 stator, 1 rotor
 "S1R2"      -- 1 stator, 2 rotor, 1 half simulated
 "S1R2_all"  -- 1 stator, 2 rotor, all simulated
 "S2R1"      -- 2 stator, 1 rotor, 1 half simulated
 "S2R1_all"  -- 2 stator, 1 rotor, all simulated

Ronald Tanner, Werner Vetter
gtisoft.com, 2023
"""
machine = {
    "name": "AFPM-10",
    "desc": "Axial flux machine",
    "afmtype": "S2R1",
    "poles": 20,
    "outer_diam": 280e-3,
    "inner_diam": 180e-3,
    "airgap": 0.003,
    "stator": {
        "num_slots": 24,
        #"mcvkey_yoke": "M330-50A",
        "afm_stator": {
            "slot_height": 0.02,
            "slot_width": 0.015,
            "slot_h1": 0.002,
            "slot_h2": 0.006,
            "slot_r1": 0.002,
            "slot_r2": 0.002,
            "yoke_height": 0.008,
            "slot_open_width": 0.003}
    },
    "magnet": {
        "nodedist": 1,
        "remanenc": 1.4,
        "afm_rotor": {
            "magn_height": 0.012,
            "rel_magn_width": 0.8,
            "yoke_height": 0,
        }
    },
    "windings": {
        "num_phases": 3,
        "num_layers": 2,
        "num_wires": 10,
        "num_par_wdgs": 1
    }
}

if __name__ == '__main__':
    import logging
    import pathlib
    import json
    import matplotlib.pyplot as plt
    from femagtools.machine.afpm import AFPM
    import femagtools.plot
    from femagtools.multiproc import Engine

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    workdir = pathlib.Path('work')
    workdir.mkdir(exist_ok=True)

    afm = AFPM(workdir, magnetizingCurves='',
               magnetMat='', condMat='')

    simulation = dict(
        angl_i_up=0.0,
        calculationMode="torq_calc",
        wind_temp=20.0,
        magn_temp=20.0,
        current=7.0711,
        speed=50.0)

    engine = Engine()
    num_slices = 3
    r = afm(engine, machine, simulation, num_slices)
    pathlib.Path('results.json').write_text(json.dumps(r))

    fig, axs = plt.subplots(nrows=2)
    femagtools.plot.torque(r['pos'], r['torque'],
                           title=f"Torque {machine['name']}", ax=axs[0])

    axs[1].set_title(f"Phase EMF (Amp {r['emf_amp']:.1f})/ V")
    axs[1].plot(r['pos'], r['emf'][0])
    axs[1].plot(r['pos'], r['emf'][1])
    axs[1].plot(r['pos'], r['emf'][2])
    axs[1].set_xlabel("Pos / °")
    axs[1].grid()
    fig.tight_layout()
    plt.show()
