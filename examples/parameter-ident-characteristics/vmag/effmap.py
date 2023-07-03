"""
  Create Efficiency Map of a VMAG PM
  Ronald Tanner
"""
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
import femagtools.plot
import femagtools.machine
import femagtools.machine.effloss

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')

    with open('dqpar.json') as fp:
        dqpars = json.load(fp)

    temp = [90, 90]
    m = femagtools.machine.create_from_eecpars(temp, dqpars)

    i1max = 250
    Udc = 420
    u1max = round(0.9*Udc/np.sqrt(2)/np.sqrt(3))
    w1, T = m.w1_imax_umax(i1max, u1max)

    nmax = 6000/60
    fig, ax = plt.subplots(figsize=(12,12))
    ax.set_title('')
    effmap = femagtools.machine.effloss.efficiency_losses_map(
        m, u1max, T, [], nmax)
    contf = femagtools.plot.efficiency_map(effmap, ax, title='')
    fig.tight_layout()
    fig.savefig('effmap.pdf')
