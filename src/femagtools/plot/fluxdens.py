""" Create airgap flux density plots

"""
import numpy as np
import matplotlib.pyplot as plt


def airgap(airgap, ax=0):
    """creates plot of flux density in airgap"""
    if ax == 0:
        ax = plt.gca()
    ax.set_title('Airgap Flux Density / T')
    ax.plot(airgap['pos'], airgap['B'],
            label='Max {:4.2f} T'.format(max(np.abs(airgap['B']))))
    ax.plot(airgap['pos'], airgap['B_fft'],
            label='Base Ampl {:4.2f} T'.format(airgap['Bamp']))
    ax.set_xlabel('Position/°')
    ax.legend()
    ax.grid(True)


def airgap_fft(airgap, bmin=1e-2, ax=0):
    """plot airgap harmonics"""
    unit = 'T'
    if ax == 0:
        ax = plt.gca()
    ax.set_title('Airgap Flux Density Harmonics / {}'.format(unit))
    ax.grid(True)
    order, fluxdens = np.array([(n, b) for n, b in zip(airgap['nue'],
                                                       airgap['B_nue']) if b > bmin]).T
    try:
        markerline1, stemlines1, _ = ax.stem(order, fluxdens, '-.', basefmt=" ")
        ax.set_xticks(order)
    except ValueError:  # empty sequence
        pass
