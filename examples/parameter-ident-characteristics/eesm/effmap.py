import logging
import json
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import femagtools.plot
import femagtools.machine
import femagtools.machine.effloss


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')

with open('eecpars.json') as fp:
    dqpars = json.load(fp)

temp = [90, 90]
m = femagtools.machine.create_from_eecpars(temp, dqpars)

i1max = 275
udc = 400
u1max = 0.9*udc/np.sqrt(3)/np.sqrt(2)
w1, Tmax = m.w1_imax_umax(i1max, u1max)
nmax = 12000/60

effmap = femagtools.machine.effloss.efficiency_losses_map(
    m, u1max, Tmax, [], nmax, num_proc=4)

fig, ax = plt.subplots(figsize=(12,12))
ax.set_title('')
contf = femagtools.plot.efficiency_map(effmap, ax, title='')
fig.tight_layout()

fig.savefig('effmap.pdf')
