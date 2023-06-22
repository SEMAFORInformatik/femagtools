import logging
import json
import matplotlib.pyplot as plt
import femagtools.plot
import femagtools.machine
import femagtools.machine.effloss


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s')


with open('impar.json') as fp:
    impars = json.load(fp)

temp = [90, 80]
T = 57
u1max = 230
nmax = 100

r = femagtools.machine.effloss.efficiency_losses_map(
    impars, u1max, T, temp, nmax, npoints=(50, 40))

fig, axs = plt.subplots(figsize=(12,10))
femagtools.plot.efficiency_map(r, ax=axs, title='')
fig.tight_layout()
fig.savefig('effmap.pdf')

# plt.show()

fig, axs = plt.subplots()
femagtools.plot.losses_map(r, ax=axs)
fig.savefig('lossmap.pdf')
# plt.show()
