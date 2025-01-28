"""Create characteristics plots

"""
import numpy as np
import matplotlib.pyplot as plt
from ..machine.utils import iqd
import logging


def mtpa(pmrel, i1max, u1max=0, title='', projection='', ax=0):
    """create a line or surface plot with torque and mtpa curve

    Args:
        pmrel: pm machine object
        i1max: maximum phase current / A
        u1max: maximum phase voltage / V (optional)
        title: plot title
        projection: projection to be used (surface if '3d')
        ax: axis to be used
    """
    nsamples = 10
    i1 = np.linspace(0, i1max, nsamples)
    iopt = np.array([pmrel.mtpa(x) for x in i1]).T

    iqmax, idmax = pmrel.iqdmax(i1max)
    iqmin, idmin = pmrel.iqdmin(i1max)

    if projection == '3d':
        nsamples = 50
    else:
        if iqmin == 0:
            iqmin = 0.1*iqmax
    id = np.linspace(idmin, idmax, nsamples)
    iq = np.linspace(iqmin, iqmax, nsamples)

    torque_iqd = np.array(
        [[pmrel.torque_iqd(x, y)
          for y in id] for x in iq])
    if projection == '3d':
        from .bch import idq_torque
        ax = idq_torque(id, iq, torque_iqd, ax)
        ax.plot(iopt[1], iopt[0], iopt[2],
                color='red', linewidth=2, label='MTPA: {0:5.0f} Nm'.format(
                    np.max(iopt[2][-1])))
    else:
        if ax == 0:
            ax = plt.gca()

        ax.set_aspect('equal')
        x, y = np.meshgrid(id, iq)
        nlevels = 6
        CST = ax.contour(x, y, torque_iqd, nlevels, colors='k')
        ax.clabel(CST, fmt='%d Nm', inline=1)
        if u1max > 0:
            n_iqd = np.array([
                [pmrel.w1_umax(u1max, iqx, idx)/np.pi/2/pmrel.p*60
                 for idx in id] for iqx in iq])
            CSN = ax.contour(x, y, n_iqd, nlevels, colors='k',
                             linestyles='dotted')  # linestyles='dashed')
            ax.clabel(CSN, fmt='%d rpm', inline=1)

        ax.set_xlabel('Id/A')
        ax.set_ylabel('Iq/A')
        ax.plot(iopt[1], iopt[0],
                color='red', linewidth=2, label='MTPA: {0:.0f} Nm'.format(
                    np.max(iopt[2][-1])))
        iqx, idx = np.array([iqd(b, i1max) for b in np.linspace(-np.pi/2, 0)]).T
        ax.plot(idx, iqx, color='blue', label='I1max: {0:.0f} A'.format(i1max))
        ax.grid()

    if title:
        ax.set_title(title)
    ax.legend(loc='upper left', framealpha=1)


def mtpv(pmrel, u1max, i1max, title='', projection='', ax=0):
    """create a line or surface plot with voltage and mtpv curve

    Args:
        pmrel: pm machine object
        u1max: maximum phase voltage / V
        i1max: maximum phase current / A
        title: plot title
        projection: projection to be used (surface if '3d')
        ax: axis to be used
    """
    w1 = pmrel.w2_imax_umax(i1max, u1max)
    nsamples = 20
    if projection == '3d':
        nsamples = 50

    iqmax, idmax = pmrel.iqdmax(i1max)
    iqmin, idmin = pmrel.iqdmin(i1max)
    id = np.linspace(idmin, idmax, nsamples)
    iq = np.linspace(iqmin, iqmax, nsamples)
    u1_iqd = np.array(
        [[np.linalg.norm(pmrel.uqd(w1, iqx, idx))/np.sqrt(2)
          for idx in id] for iqx in iq])
    u1 = np.mean(u1_iqd)
    imtpv = np.array([pmrel.mtpv(wx, u1, i1max)
                      for wx in np.linspace(w1, 20*w1, nsamples)]).T

    if projection == '3d':
        from .bch import idq_torque
        torque_iqd = np.array(
            [[pmrel.torque_iqd(x, y)
              for y in id] for x in iq])
        ax = idq_torque(id, iq, torque_iqd, ax)
        ax.plot(imtpv[1], imtpv[0], imtpv[2],
                color='red', linewidth=2)
    else:
        if ax == 0:
            ax = plt.gca()
        ax.set_aspect('equal')
        x, y = np.meshgrid(id, iq)
        CS = ax.contour(x, y, u1_iqd, 4, colors='b')  # linestyles='dashed')
        ax.clabel(CS, fmt='%d', inline=1)

        ax.plot(imtpv[1], imtpv[0],
                color='red', linewidth=2,
                label='MTPV: {0:5.0f} Nm'.format(np.max(imtpv[2])))
        # beta = np.arctan2(imtpv[1][0], imtpv[0][0])
        # b = np.linspace(beta, 0)
        # ax.plot(np.sqrt(2)*i1max*np.sin(b), np.sqrt(2)*i1max*np.cos(b), 'r-')

        ax.grid()
        ax.legend()
    ax.set_xlabel('Id/A')
    ax.set_ylabel('Iq/A')
    if title:
        ax.set_title(title)


def characteristics(char, title=''):
    """show 2x2 line plot as speed characteristics: (torque, pmech), (u1, cosphi)
    (i1, id, iq, beta), (fe, cu, fric loss, eta)

    Args:
        char: result from characteristics calculation (im, sm, pm)
        title: figure title

    Returns:
        fig object
    """
    fig, axs = plt.subplots(2, 2, figsize=(10, 8),
                            layout='tight',
                            sharex=True)
    if title:
        fig.suptitle(title)

    n = np.array(char['n'])*60
    punit = 'kW'
    k = 1e-3
    if max(char['pmech']) > 1e6:
        punit = 'MW'
        k = 1e-6
    pmech = np.array(char['pmech'])*k
    tunit = 'Nm'
    if max(char['T']) > 1e3:
        tunit = 'kNm'
        tq = np.array(char['T'])*1e-3
    else:
        tq = np.array(char['T'])

    axs[0, 0].plot(n, tq, 'C0-', label='Torque')
    axs[0, 0].set_ylabel(f"Torque / {tunit}")
    axs[0, 0].grid()
    axs[0, 0].legend(loc='center left')
    ax1 = axs[0, 0].twinx()
    ax1.plot(n, pmech, 'C1-', label='P mech')
    ax1.set_ylabel(f"Power / {punit}")
    ax1.legend(loc='lower center')

    axs[0, 1].plot(n[1:], np.array(char['u1'][1:]), 'C0-', label='Voltage')
    axs[0, 1].set_ylabel("Phase Voltage / V",)
    axs[0, 1].grid()
    axs[0, 1].legend(loc='center left')
    ax2 = axs[0, 1].twinx()
    ax2.plot(n[1:], char['cosphi'][1:], 'C1-', label='Cos Phi')
    ax2.set_ylabel("Cos Phi")
    ax2.legend(loc='lower right')

    if 'id' in char:
        axs[1, 0].plot(n, np.array(char['id']), label='Id')
    if 'iq' in char:
        axs[1, 0].plot(n, np.array(char['iq']), label='Iq')
    axs[1, 0].plot(n, np.array(char['i1']), label='I1')
    axs[1, 0].set_xlabel("Speed / rpm")
    axs[1, 0].set_ylabel("Current / A")
    axs[1, 0].legend(loc='center left')
    if 'iex' in char:
        ax3 = axs[1, 0].twinx()
        ax3.plot(n, char['iex'], 'C3-', label='Iexc')
        ax3.set_ylabel("Iexc / A")
        ax3.legend(loc='center right')
    elif 'beta' in char:
        ax3 = axs[1, 0].twinx()
        ax3.plot(n, char['beta'], 'C3-', label='Beta')
        ax3.set_ylabel("Beta / °")
        ax3.legend(loc='center right')
    axs[1, 0].grid()
    try:
        plfe = np.array(char['plfe'])*1e-3
    except KeyError:
        plfe = np.array(char['plfe1'])*1e-3
    try:
        plcu = np.array(char['plcu'])*1e-3
    except KeyError:
        plcu = np.array(char['plcu1'])*1e-3
    pl = np.array(char['losses'])*1e-3
    axs[1, 1].plot(n, plcu, 'C0-', label='Cu Losses')
    axs[1, 1].plot(n, plfe, 'C1-', label='Fe Losses')
    try:
        if char['plfw'] and char['plfw'][-1] > 0:
            plfw = np.array(char['plfw'])*1e-3
            axs[1, 1].plot(n, plfw, 'C2-', label='Friction + Windage')
    except KeyError:
            pass
    axs[1, 1].set_ylabel("Losses / kW")
    axs[1, 1].legend(loc='center left')
    axs[1, 1].grid()
    axs[1, 1].set_xlabel("Speed / rpm")
    ax4 = axs[1, 1].twinx()
    ax4.plot(n[1:-1], char['eta'][1:-1], 'C3-', label="Eta")
    ax4.legend(loc='upper center')
    ax4.set_ylabel("Efficiency")

    return fig


def _get_nT_boundary(n, T):
    """utility function to extract boundary from n,T lists"""
    n0 = n[0]
    t0 = T[0]
    #      braking, driving
    bnd = [[(n0, t0)], []]
    for nx, tx in zip(n, T):
        if tx < t0:
            bnd[1].append((n0, t0))
            bnd[0].append((nx, tx))
            n0 = nx
        t0 = tx
    bnd[1].append((nx, tx))
    return np.array(bnd[0] + bnd[1][::-1])

def _normalize10(v, **kwargs):
    """utility function to normalize the input-array using the nearest (ceiling) power of 10.

    Args:
        v: array to normalize

    Returns:
        normalized array
        normalisation factor (power of 10)
    """
    n_type = kwargs.get('n_type', 'abs') # 'norm')
    r_type = kwargs.get('r_type', 'round') # 'floor')
    if n_type == 'norm':
        norm = np.log10(np.linalg.norm(v))
    elif n_type == 'abs':
        norm = np.log10(np.max(np.abs(v)))
    else:
        raise AttributeError('Unknown norm-type. Allowed values are: "norm", "abs".')

    if r_type == 'floor':
        norm = 10 ** np.floor(norm)
    elif r_type == 'ceil':
        norm = 10 ** np.ceil(norm)
    elif r_type == 'round':
        norm = 10 ** np.round(norm)
    else:
        raise AttributeError('Unknown rounding type. Allowed values are: "floor", "ceil", "round".')
    # norm = 10**(np.ceil(np.log10(np.max(np.abs(v)))))
    # norm = 10 ** (np.ceil(np.log10(np.linalg.norm(v))))
    # norm = 10 ** (np.floor(np.log10(np.max(np.abs(v)))))
    # norm = 10 ** (np.floor(np.log10(np.linalg.norm(v))))
    return v / norm, norm


def _plot_contour(speed, torque, z, ax, title='', levels=[],
                  clabel=True, cmap='YlOrRd', cbar=False, **kwargs):
    """utility function for contour plot of speed, torque, z values

    Args:
        levels: (list of floats)
        clabel: contour labels if True
        cmap: colour map
        cbar: (bool) create color bar if True (default False)

    Note: the x and y axes are scaled

    Returns:
        tricontourf, xscale, yscale
    """
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    x = 60*np.asarray(speed)
    y = np.asarray(torque)
    tunit = 'Nm'
    if max(np.abs(y)) > 10e3:
        y *=1e-3
        tunit = 'kNm'
    x, xscale = _normalize10(x, **kwargs)
    y, yscale = _normalize10(y, **kwargs)

    if not levels:
        if max(z) <= 1:
            if max(z) > 0.96:
                levels = [0.5, 0.75, 0.8, 0.84,
                          0.89, 0.92, 0.94, 0.96, 0.97]
            else:
                levels = [0.25, 0.5, 0.75, 0.8, 0.84,
                          0.88, 0.9, 0.92, 0.94, 0.96]

            if max(z) > levels[-1]:
                levels.append(np.ceil(max(z)*100)/100)
        else:
            levels = 14
    #
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')

    cont = ax.tricontour(x, y, z,
                         linewidths=0.4, levels=levels,
                         colors='k')
    contf = ax.tricontourf(x, y, z,
                           levels=levels, cmap=cmap)

    #ax.plot(x, y, "k.", ms=3)
    if clabel:
        ax.clabel(cont, inline=True, colors='k', fontsize=8, inline_spacing=0)

    clippath = Path(_get_nT_boundary(x, y))
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    try:
        for c in cont.collections:
            c.set_clip_path(patch)
        for c in contf.collections:
            c.set_clip_path(patch)
    except AttributeError:  # matplotlib >= 3.10
        cont.set_clip_path(patch)
        contf.set_clip_path(patch)

    if xscale > 1:
        def format_fn(tick_val, tick_pos):
            return round(xscale*tick_val)
        ax.xaxis.set_major_formatter(format_fn)
    if yscale > 1:
        def format_fn(tick_val, tick_pos):
            return round(yscale*tick_val)
        ax.yaxis.set_major_formatter(format_fn)

    ax.set_ylabel(f'Torque / {tunit}')
    ax.set_xlabel('Speed / rpm')
    ax.set_title(title)
    if cbar:
        cfig = ax.get_figure()
        cfig.colorbar(contf, ax=ax,
                      orientation='vertical')
    return contf, xscale, yscale


def efficiency_map(rmap, ax=0, title='', clabel=True,
                   cmap='YlOrRd', levels=None, cbar=False):
    if ax == 0:
        fig, ax = plt.subplots(figsize=(12, 12))
    return _plot_contour(rmap['n'], rmap['T'], rmap['eta'], ax,
                        title=title, clabel=clabel, cmap=cmap,
                        levels=levels, cbar=cbar)


def losses_map(rmap, ax=0, title='Losses Map / kW', clabel=True,
               cmap='YlOrRd', cbar=False, key='losses'):
    """
    plot losses map
    Args:
    rmap: (dict) result of efficiency_losses_map
    key: (str) type of losses: 'plfe1', 'plfe2', 'plmag', 'plcu1', 'plcu2', 'plfric', 'losses';
    """

    if ax == 0:
        fig, ax = plt.subplots(figsize=(12, 12))
    return _plot_contour(rmap['n'], rmap['T'], np.asarray(rmap[key])/1e3, ax,
                         title=title, levels=14, clabel=clabel,
                         cmap=cmap, cbar=cbar)
