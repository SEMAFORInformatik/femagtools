import logging
import json
import matplotlib.pyplot as plt
import numpy as np

import femagtools.plot
import femagtools.machine
import femagtools.machine.effloss


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

with open('eecpars.json') as fp:
    dqpars = json.load(fp)

temp = [90, 90]
m = femagtools.machine.create_from_eecpars(temp, dqpars)

T = 240
udc = 400
u1max = 0.9*udc/np.sqrt(3)/np.sqrt(2)
nmax = 12000/60
fig, ax = plt.subplots(figsize=(12,12))
ax.set_title('')
effmap = femagtools.machine.effloss.efficiency_losses_map(
    m, u1max, T, [], nmax)
contf = femagtools.plot.efficiency_map(effmap, ax, title='')
fig.tight_layout()
##fig.colorbar(contf)
fig.savefig('effmap.pdf')
