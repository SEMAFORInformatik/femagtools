"""
    femagtools.plot.bch
    ~~~~~~~~~~~~~~~~~~~

    Creating bch plots


"""
import numpy as np
import scipy.interpolate as ip
import logging

try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from mpl_toolkits.mplot3d import Axes3D
    matplotlibversion = matplotlib.__version__
except ImportError:   # ModuleNotFoundError:
    matplotlibversion = 0

logger = logging.getLogger("femagtools.plot.bch")


def find_peaks_and_valleys(t, y):
    """ return peaks and valleys of y with maximum amplitude
    """
    peaks = (np.diff(np.sign(np.diff(y))) < 0).nonzero()[0] + 1
    if len(peaks > 0):
        ip = np.argmax(y[peaks])
        pv = {'yp': y[peaks][ip], 'tp': t[peaks][ip]}
    else:
        pv = {'yp': [], 'tp': []}
    valleys = (np.diff(np.sign(np.diff(y))) > 0).nonzero()[0] + 1
    if len(valleys > 0):
        iv = np.argmin(y[valleys])
        pv.update({'yv': y[valleys][iv], 'tv': t[valleys][iv]})
    else:
        pv.update({'yv': [], 'tv': []})
    pv.update({'peaks': y[peaks], 'valleys': y[valleys],
               'tpeaks': t[peaks], 'tvalleys': t[valleys]})
    return pv


def _create_3d_axis():
    """creates a subplot with 3d projection if one does not already exist"""
    from matplotlib.projections import get_projection_class
    from matplotlib import _pylab_helpers

    create_axis = True
    if _pylab_helpers.Gcf.get_active() is not None:
        if isinstance(plt.gca(), get_projection_class('3d')):
            create_axis = False
    if create_axis:
        plt.figure()
        plt.subplot(111, projection='3d')


def _plot_surface(ax, x, y, z, labels, azim=None):
    """helper function for surface plots"""
    # ax.tick_params(axis='both', which='major', pad=-3)
    assert np.size(x) > 1 and np.size(y) > 1 and np.size(z) > 1
    if azim is not None:
        ax.azim = azim
    X, Y = np.meshgrid(x, y)
    Z = np.ma.masked_invalid(z)
    ax.plot_surface(X, Y, Z,
                    rstride=1, cstride=1,
                    cmap=cm.viridis, alpha=0.85,
                    vmin=np.nanmin(z), vmax=np.nanmax(z),
                    linewidth=0, antialiased=True)
#                    edgecolor=(0, 0, 0, 0))
    # ax.set_xticks(xticks)
    # ax.set_yticks(yticks)
    # ax.set_zticks(zticks)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(labels[2])

    # plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)


def torque(pos, torque, title='', ax=0):
    """creates plot from torque vs position"""
    k = 20
    alpha = np.linspace(pos[0], pos[-1],
                        k*len(torque))
    f = ip.interp1d(pos, torque, kind='quadratic')
    unit = 'Nm'
    scale = 1
    if np.min(torque) < -9.9e3 or np.max(torque) > 9.9e3:
        scale = 1e-3
        unit = 'kNm'
    if ax == 0:
        ax = plt.gca()
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Torque / {}'.format(unit))
    ax.grid(True)
    ax.plot(pos, [scale*t for t in torque], 'go')
    ax.plot(alpha, scale*f(alpha))
    if np.min(torque) > 0 and np.max(torque) > 0:
        ax.set_ylim(bottom=0)
    elif np.min(torque) < 0 and np.max(torque) < 0:
        ax.set_ylim(top=0)


def torque_fft(order, torque, ax=0):
    """plot torque harmonics"""
    unit = 'Nm'
    scale = 1
    if np.min(torque) < -9.9e3 or np.max(torque) > 9.9e3:
        scale = 1e-3
        unit = 'kNm'
    if ax == 0:
        ax = plt.gca()
    ax.set_title('Torque Harmonics / {}'.format(unit))
    ax.grid(True)

    try:
        bw = 2.5E-2*max(order)
        ax.bar(order, [scale*t for t in torque], width=bw, align='center')
        ax.set_xlim(left=-bw/2)
    except ValueError:  # empty sequence
        pass


def force(title, pos, force, xlabel='', ax=0):
    """plot force vs position"""
    unit = 'N'
    scale = 1
    if min(force) < -9.9e3 or max(force) > 9.9e3:
        scale = 1e-3
        unit = 'kN'
    if ax == 0:
        ax = plt.gca()
    ax.set_title('{} / {}'.format(title, unit))
    ax.grid(True)
    ax.plot(pos, [scale*f for f in force])
    if xlabel:
        ax.set_xlabel(xlabel)
    if min(force) > 0:
        ax.set_ylim(bottom=0)


def force_fft(order, force, ax=0):
    """plot force harmonics"""
    unit = 'N'
    scale = 1
    if min(force) < -9.9e3 or max(force) > 9.9e3:
        scale = 1e-3
        unit = 'kN'
    if ax == 0:
        ax = plt.gca()
    ax.set_title('Force Harmonics / {}'.format(unit))
    ax.grid(True)
    try:
        bw = 2.5E-2*max(order)
        ax.bar(order, [scale*t for t in force], width=bw, align='center')
        ax.set_xlim(left=-bw/2)
    except ValueError:  # empty sequence
        pass


def fluxdens_surface(fdens, ax=0):
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    x = [p for p in fdens.positions[0]['X']]
    y = [p['position'] for p in fdens.positions]
    z = np.array([p['B_N']
                  for p in fdens.positions])
    _plot_surface(ax, x, y, z,
                  (u'Rotor pos/°', u'Pos/°', u'B N / T'))


def winding_flux(pos, flux, ax=0):
    """plot flux vs position"""
    if ax == 0:
        ax = plt.gca()
    ax.set_title('Winding Flux / Vs')
    ax.grid(True)
    for p, f in zip(pos, flux):
        ax.plot(p, f)


def winding_current(pos, current, ax=0):
    """plot winding currents"""
    if ax == 0:
        ax = plt.gca()
    ax.set_title('Winding Currents / A')
    ax.grid(True)
    for p, i in zip(pos, current):
        ax.plot(p, i)


def voltage(title, pos, voltage, ax=0):
    """plot voltage vs. position"""
    if ax == 0:
        ax = plt.gca()
    ax.set_title('{} / V'.format(title))
    ax.grid(True)
    ax.plot(pos, voltage)


def voltage_fft(title, order, voltage, ax=0):
    """plot FFT harmonics of voltage"""
    if ax == 0:
        ax = plt.gca()
    ax.set_title('{} / V'.format(title))
    ax.grid(True)
    if max(order) < 5:
        order += [5]
        voltage += [0]
    try:
        bw = 2.5E-2*max(order)
        ax.bar(order, voltage, width=bw, align='center')
    except ValueError:  # empty sequence
        pass


def __get_linearForce_title_keys(lf):
    if 'force_r' in lf:
        return ['Force r', 'Force z'], ['force_r', 'force_z']
    return ['Force x', 'Force y'], ['force_x', 'force_y']


def pmrelsim(bch, title=''):
    """creates a plot of a PM/Rel motor simulation"""
    cols = 2
    rows = 4
    if len(bch.flux['1']) > 1:
        rows += 1
    htitle = 1.5 if title else 0
    fig, ax = plt.subplots(nrows=rows, ncols=cols,
                           figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)

    row = 1
    plt.subplot(rows, cols, row)
    if bch.torque:
        torque(bch.torque[-1]['angle'], bch.torque[-1]['torque'])
        plt.subplot(rows, cols, row+1)
        tq = list(bch.torque_fft[-1]['torque'])
        order = list(bch.torque_fft[-1]['order'])
        if order and max(order) < 5:
            order += [15]
            tq += [0]
        torque_fft(order, tq)
        plt.subplot(rows, cols, row+2)
        force('Force Fx',
              bch.torque[-1]['angle'], bch.torque[-1]['force_x'])
        plt.subplot(rows, cols, row+3)
        force('Force Fy',
              bch.torque[-1]['angle'], bch.torque[-1]['force_y'])
        row += 3
    elif bch.linearForce:
        title, keys = __get_linearForce_title_keys(bch.linearForce[-1])
        force(title[0], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[0]], 'Displt. / mm')
        plt.subplot(rows, cols, row+1)
        force_fft(bch.linearForce_fft[-2]['order'],
                  bch.linearForce_fft[-2]['force'])
        plt.subplot(rows, cols, row+2)
        force(title[1], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[1]], 'Displt. / mm')
        plt.subplot(rows, cols, row+3)
        force_fft(bch.linearForce_fft[-1]['order'],
                  bch.linearForce_fft[-1]['force'])
        row += 3

    plt.subplot(rows, cols, row+1)
    flux = [bch.flux[k][-1] for k in bch.flux]
    pos = [f['displ'] for f in flux]
    winding_flux(pos,
                 [f['flux_k'] for f in flux])
    plt.subplot(rows, cols, row+2)
    winding_current(pos,
                    [f['current_k'] for f in flux])
    plt.subplot(rows, cols, row+3)
    voltage('Internal Voltage',
            bch.flux['1'][-1]['displ'],
            bch.flux['1'][-1]['voltage_dpsi'])
    plt.subplot(rows, cols, row+4)
    try:
        voltage_fft('Internal Voltage Harmonics',
                    bch.flux_fft['1'][-1]['order'],
                    bch.flux_fft['1'][-1]['voltage'])
    except:
        pass
    if len(bch.flux['1']) > 1:
        plt.subplot(rows, cols, row+5)
        voltage('No Load Voltage',
                bch.flux['1'][0]['displ'],
                bch.flux['1'][0]['voltage_dpsi'])
        plt.subplot(rows, cols, row+6)
        try:
            voltage_fft('No Load Voltage Harmonics',
                        bch.flux_fft['1'][0]['order'],
                        bch.flux_fft['1'][0]['voltage'])
        except:
            pass
    fig.tight_layout(h_pad=3.5)
    if title:
        fig.subplots_adjust(top=0.92)


def multcal(bch, title=''):
    """creates a plot of a MULT CAL simulation"""
    cols = 2
    rows = 4
    htitle = 1.5 if title else 0
    fig, ax = plt.subplots(nrows=rows, ncols=cols,
                           figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)

    row = 1
    plt.subplot(rows, cols, row)
    if bch.torque:
        torque(bch.torque[-1]['angle'], bch.torque[-1]['torque'])
        plt.subplot(rows, cols, row+1)
        tq = list(bch.torque_fft[-1]['torque'])
        order = list(bch.torque_fft[-1]['order'])
        if order and max(order) < 5:
            order += [15]
            tq += [0]
        torque_fft(order, tq)
        plt.subplot(rows, cols, row+2)
        force('Force Fx',
              bch.torque[-1]['angle'], bch.torque[-1]['force_x'])
        plt.subplot(rows, cols, row+3)
        force('Force Fy',
              bch.torque[-1]['angle'], bch.torque[-1]['force_y'])
        row += 3
    elif bch.linearForce:
        title, keys = __get_linearForce_title_keys(bch.linearForce[-1])
        force(title[0], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[0]], 'Displt. / mm')
        plt.subplot(rows, cols, row+1)
        force_fft(bch.linearForce_fft[-2]['order'],
                  bch.linearForce_fft[-2]['force'])
        plt.subplot(rows, cols, row+2)
        force(title[1], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[1]], 'Displt. / mm')
        plt.subplot(rows, cols, row+3)
        force_fft(bch.linearForce_fft[-1]['order'],
                  bch.linearForce_fft[-1]['force'])
        row += 3

    plt.subplot(rows, cols, row+1)
    flux = [bch.flux[k][-1] for k in bch.flux]
    pos = [f['displ'] for f in flux]
    winding_flux(pos,
                 [f['flux_k'] for f in flux])
    plt.subplot(rows, cols, row+2)
    winding_current(pos,
                    [f['current_k'] for f in flux])
    plt.subplot(rows, cols, row+3)
    voltage('Internal Voltage',
            bch.flux['1'][-1]['displ'],
            bch.flux['1'][-1]['voltage_dpsi'])
    plt.subplot(rows, cols, row+4)
    try:
        voltage_fft('Internal Voltage Harmonics',
                    bch.flux_fft['1'][-1]['order'],
                    bch.flux_fft['1'][-1]['voltage'])
    except:
        pass
    if len(bch.flux['1']) > 1:
        plt.subplot(rows, cols, row+5)
        voltage('No Load Voltage',
                bch.flux['1'][0]['displ'],
                bch.flux['1'][0]['voltage_dpsi'])
        plt.subplot(rows, cols, row+6)
        try:
            voltage_fft('No Load Voltage Harmonics',
                        bch.flux_fft['1'][0]['order'],
                        bch.flux_fft['1'][0]['voltage'])
        except:
            pass

    fig.tight_layout(h_pad=3.5)
    if title:
        fig.subplots_adjust(top=0.92)


def fasttorque(bch, title=''):
    """creates a plot of a Fast Torque simulation"""
    cols = 2
    rows = 4
    if len(bch.flux['1']) > 1:
        rows += 1
    htitle = 1.5 if title else 0
    fig, ax = plt.subplots(nrows=rows, ncols=cols,
                           figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)

    row = 1
    plt.subplot(rows, cols, row)
    if bch.torque:
        torque(bch.torque[-1]['angle'], bch.torque[-1]['torque'])
        plt.subplot(rows, cols, row+1)
        torque_fft(bch.torque_fft[-1]['order'], bch.torque_fft[-1]['torque'])
        plt.subplot(rows, cols, row+2)
        force('Force Fx',
              bch.torque[-1]['angle'], bch.torque[-1]['force_x'])
        plt.subplot(rows, cols, row+3)
        force('Force Fy',
              bch.torque[-1]['angle'], bch.torque[-1]['force_y'])
        row += 3
    elif bch.linearForce:
        title, keys = __get_linearForce_title_keys(bch.linearForce[-1])
        force(title[0], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[0]], 'Displt. / mm')
        plt.subplot(rows, cols, row+1)
        force_fft(bch.linearForce_fft[-2]['order'],
                  bch.linearForce_fft[-2]['force'])
        plt.subplot(rows, cols, row+2)
        force(title[1], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[1]], 'Displt. / mm')
        plt.subplot(rows, cols, row+3)
        force_fft(bch.linearForce_fft[-1]['order'],
                  bch.linearForce_fft[-1]['force'])
        row += 3

    plt.subplot(rows, cols, row+1)
    flux = [bch.flux[k][-1] for k in bch.flux]
    pos = [f['displ'] for f in flux]
    winding_flux(pos, [f['flux_k'] for f in flux])
    plt.subplot(rows, cols, row+2)
    winding_current(pos, [f['current_k'] for f in flux])
    plt.subplot(rows, cols, row+3)
    voltage('Internal Voltage',
            bch.flux['1'][-1]['displ'],
            bch.flux['1'][-1]['voltage_dpsi'])
    plt.subplot(rows, cols, row+4)
    try:
        voltage_fft('Internal Voltage Harmonics',
                    bch.flux_fft['1'][-1]['order'],
                    bch.flux_fft['1'][-1]['voltage'])
    except:
        pass
    if len(bch.flux['1']) > 1:
        plt.subplot(rows, cols, row+5)
        voltage('No Load Voltage',
                bch.flux['1'][0]['displ'],
                bch.flux['1'][0]['voltage_dpsi'])
        plt.subplot(rows, cols, row+6)
        try:
            voltage_fft('No Load Voltage Harmonics',
                        bch.flux_fft['1'][0]['order'],
                        bch.flux_fft['1'][0]['voltage'])
        except:
            pass
    fig.tight_layout(h_pad=3.5)
    if title:
        fig.subplots_adjust(top=0.92)


def cogging(bch, title=''):
    """creates a cogging plot"""
    cols = 2
    rows = 3

    htitle = 1.5 if title else 0
    fig, ax = plt.subplots(nrows=rows, ncols=cols,
                           figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)

    row = 1
    plt.subplot(rows, cols, row)
    if bch.torque:
        torque(bch.torque[0]['angle'], bch.torque[0]['torque'])
        plt.subplot(rows, cols, row+1)
        if bch.torque_fft:
            torque_fft(bch.torque_fft[0]['order'], bch.torque_fft[0]['torque'])
        plt.subplot(rows, cols, row+2)
        force('Force Fx',
              bch.torque[0]['angle'], bch.torque[0]['force_x'])
        plt.subplot(rows, cols, row+3)
        force('Force Fy',
              bch.torque[0]['angle'], bch.torque[0]['force_y'])
        row += 3
    elif bch.linearForce:
        title, keys = __get_linearForce_title_keys(bch.linearForce[-1])
        force(title[0], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[0]], 'Displt. / mm')
        plt.subplot(rows, cols, row+1)
        force_fft(bch.linearForce_fft[-2]['order'],
                  bch.linearForce_fft[-2]['force'])
        plt.subplot(rows, cols, row+2)
        force(title[1], bch.linearForce[-1]['displ'],
              bch.linearForce[-1][keys[1]], 'Displt. / mm')
        plt.subplot(rows, cols, row+3)
        force_fft(bch.linearForce_fft[-1]['order'],
                  bch.linearForce_fft[-1]['force'])
        row += 3

    plt.subplot(rows, cols, row+1)
    voltage('Voltage',
            bch.flux['1'][0]['displ'],
            bch.flux['1'][0]['voltage_dpsi'])
    plt.subplot(rows, cols, row+2)
    voltage_fft('Voltage Harmonics',
                bch.flux_fft['1'][0]['order'],
                bch.flux_fft['1'][0]['voltage'])

    fig.tight_layout(h_pad=2)
    if title:
        fig.subplots_adjust(top=0.92)


def demagnetization(demag, ax=0):
    """plot rel. remanence vs. current"""
    if ax == 0:
        ax = plt.gca()
    scale = 1
    unit = 'A'
    if np.max(demag['i1']) > 25e3:
        scale = 1e-3
        unit = 'kA'
    i1 = [scale*x for x in demag['i1']]
    ax.plot(i1, demag['rr'], 'o', color='C0')
    ax.plot(i1, demag['rr'], color='C0')
    rrmin = 0.6
    if demag.get('i1c', 0):
        Icrit = scale*demag['i1c']
        Hk = demag['Hk']
        Tmag = demag['Tmag']
        di = 0.05*np.max(i1)
        rrmin = min(0.6, np.min(demag['rr']))
        ax.plot([Icrit, Icrit], [rrmin, 1], 'k--')
        ax.annotate(
            f'Icrit = {Icrit:.1f}{unit}\nHk = {Hk:.1f} kA/m\nTmag={Tmag:.1f} °C',
            xy=(Icrit, rrmin),
            xytext=(Icrit-di, rrmin+0.1*(1-rrmin)), ha='right',
            bbox={'facecolor': 'white',
                  'edgecolor': 'white'})
    ax.set_ylim([rrmin, 1.01])
    ax.set_ylabel('Rel. Remanence')
    ax.set_xlabel(f'Phase Current / {unit}')
    ax.grid()


def transientsc_currents(scData, ax=0, title='', set_xlabel=True):
    """plot transient shortcircuit currents vs time"""
    if ax == 0:
        ax = plt.gca()
    ax.grid(True)
    if title:
        ax.set_title(title)
    istat = np.array([scData[i]
                      for i in ('ia', 'ib', 'ic')])
    pv = [find_peaks_and_valleys(
        np.array(scData['time']), i1)
          for i1 in istat]
    try:
        ipvmax = np.argmax(
            [y['yp'] if np.abs(y['yp']) > np.abs(y['yv']) else y['yv']
             for y in pv if y['yp']])
        imax = pv[ipvmax]['yp'] if np.abs(pv[ipvmax]['yp']) > np.abs(pv[ipvmax]['yv']) else pv[ipvmax]['yv']
        iac = [pv[ipvmax]['tpeaks'][-1], pv[ipvmax]['peaks'][-1]]
    except KeyError:
        pass
    if np.max(istat) > 4000:
        istat *= 1e-3
        try:
            imax *= 1e-3
            iac[1] *= 1e-3
        except NameError:
            pass
        ax.set_ylabel('Currents / kA')
    else:
        ax.set_ylabel('Currents / A')

    for i, iph in zip(('ia', 'ib', 'ic'), istat):
        ax.plot(scData['time'], iph, label=i)
    try:
        ax.plot([pv[ipvmax]['tp']], [imax], '.')
        ax.plot([iac[0]], [iac[1]], '.')
        dtx = (scData['time'][-1]-scData['time'][0])/75
        dy = imax/25
        ax.annotate(f'Imax = {imax:.1f}',
                    xy=(pv[ipvmax]['tp'], imax),
                    xytext=(pv[ipvmax]['tp']+dtx, imax-dy))
        dy = iac[1]/25
        ax.annotate(f'I = {iac[1]:.1f}',
                    xy=iac,
                    xytext=(iac[0]+dtx, iac[1]-dy))
    except NameError:
        pass
    if set_xlabel:
        ax.set_xlabel('Time / s')
    ax.legend()


def transientsc_torque(scData, ax=0, title='', set_xlabel=True):
    """plot transient shortcircuit torque vs time"""
    if ax == 0:
        ax = plt.gca()
    if title:
        ax.set_title(title)
    pv = find_peaks_and_valleys(
        np.array(scData['time']), np.array(scData['torque']))
    try:
        tqmax = pv['yp'] if np.abs(pv['yp']) > np.abs(pv['yv']) else pv['yv']
        tp = pv['tp'] if np.abs(pv['yp']) > np.abs(pv['yv']) else pv['tv']
        tc = [pv['tpeaks'][-1], pv['peaks'][-1]]
    except (KeyError, ValueError):
        pass
    torque = np.array(scData['torque'])
    if np.max(torque) > 4000:
        torque *= 1e-3
        try:
            tqmax *= 1e-3
            tc[1] *= 1e-3
        except NameError:
            pass
        ax.set_ylabel('Torque / kNm')
    else:
        ax.set_ylabel('Torque / Nm')

    ax.grid(True)
    ax.plot(scData['time'], torque)
    try:
        ax.plot([tp], [tqmax], '.')
        ax.plot([tc[0]], [tc[1]], '.')
        dtx = (scData['time'][-1]-scData['time'][0])/75
        dy = tqmax/25
        ax.annotate(f'Tmax = {tqmax:.1f}',
                    xy=(tp, tqmax),
                    xytext=(tp+dtx, tqmax-dy))
        dy = tc[1]/25
        ax.annotate(f'T = {tc[1]:.1f}',
                    xy=tc,
                    xytext=(tc[0]+dtx, tc[1]))
    except NameError:
        pass
    if set_xlabel:
        ax.set_xlabel('Time / s')

def transientsc(bch, title=''):
    """creates a transient short circuit plot"""
    try:
        scData = bch.scData
    except AttributeError:
        scData = bch
    cols = 1
    rows = 2
    htitle = 1.5 if title else 0
    fig, axs = plt.subplots(nrows=rows, ncols=cols, sharex=True,
                            figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)

    transientsc_currents(scData, axs[0], set_xlabel=False)
    transientsc_torque(scData, axs[1])
    fig.tight_layout(h_pad=2)
    if title:
        fig.subplots_adjust(top=0.92)
    return fig

def transientsc_demag(demag, magnet=0, title='', ax=0):
    """creates a demag plot of a transient short circuit
    Args:
      demag: list of dicts with 'displ', 'H_av', 'H_max', 'lim_hc'
      magnet dict with 'Tmag'
    """
    if ax == 0:
        ax = plt.gca()
    if type(d) == list:
        pos = [d['displ'] for d in demag if 'displ' in d]
        hmax = [-d['H_max'] for d in demag if 'H_max' in d]
        havg = [-d['H_av'] for d in demag if 'H_av' in d]
        hclim = [-d['lim_hc'] for d in demag if 'lim_hc' in d][0]
    else:
        pos = demag['displ']
        hmax = demag['H_max']
        havg = demag['H_av']
        hclim = demag['Hk']
    ax.set_title('Transient Short Circuit Demagnetization [kA/m]')
    ax.plot(pos, hmax,
            label='H Max {:4.2f} kA/m'.format(max(hmax)))
    ax.plot(pos, havg,
            label='H Avg {:4.2f} kA/m'.format(max(havg)))
    if len(hclim) > 1:
        ax.plot([pos[0], pos[-1]], [hclim,hclim], color='C3',
                linestyle='dashed',
                label='Hc {:4.2f} kA/m'.format(hclim))
    if 'Tmag' in demag:
        Tmag = demag['Tmag']
    elif magnet:
        Tmag = magnet['Tmag']
    else:
        Tmag = ''
    ax.set_xlabel('Rotor Position / °')
    ax.grid(True)
    if Tmag:
        ax.legend(title=f"Magnet Temperature {Tmag}°C")
    else:
        ax.legend()

def i1beta_torque(i1, beta, torque, title='', ax=0):
    """creates a surface plot of torque vs i1, beta"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    azim = 210
    if 0 < np.mean(beta) or -90 > np.mean(beta):
        azim = -60
    unit = 'Nm'
    scale = 1
    if np.min(torque) < -9.9e3 or np.max(torque) > 9.9e3:
        scale = 1e-3
        unit = 'kNm'
    if title:
        _plot_surface(ax, i1, beta, scale*np.asarray(torque),
                      (u'I1/A', u'Beta/°', title),
                      azim=azim)
    else:
        _plot_surface(ax, i1, beta, scale*np.asarray(torque),
                      (u'I1/A', u'Beta/°', u'Torque/{}'.format(unit)),
                      azim=azim)


def i1beta_ld(i1, beta, ld, ax=0):
    """creates a surface plot of ld vs i1, beta"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, i1, beta, np.asarray(ld)*1e3,
                  (u'I1/A', u'Beta/°', u'Ld/mH'),
                  azim=60)


def i1beta_lq(i1, beta, lq, ax=0):
    """creates a surface plot of ld vs i1, beta"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    azim = 60
    if 0 < np.mean(beta) or -90 > np.mean(beta):
        azim = -120
    _plot_surface(ax, i1, beta, np.asarray(lq)*1e3,
                  (u'I1/A', u'Beta/°', u'Lq/mH'),
                  azim=azim)


def i1beta_psim(i1, beta, psim, ax=0):
    """creates a surface plot of psim vs i1, beta"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, i1, beta, psim,
                  (u'I1/A', u'Beta/°', u'Psi m/Vs'),
                  azim=60)


def i1beta_up(i1, beta, up, ax=0):
    """creates a surface plot of up vs i1, beta"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, i1, beta, up,
                  (u'I1/A', u'Beta/°', u'Up/V'),
                  azim=60)


def i1beta_psid(i1, beta, psid, ax=0):
    """creates a surface plot of psid vs i1, beta"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    azim = -60
    if 0 < np.mean(beta) or -90 > np.mean(beta):
        azim = 60
    _plot_surface(ax, i1, beta, psid,
                  (u'I1/A', u'Beta/°', u'Psi d/Vs'),
                  azim=azim)


def i1beta_psiq(i1, beta, psiq, ax=0):
    """creates a surface plot of psiq vs i1, beta"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    azim = 210
    if 0 < np.mean(beta) or -90 > np.mean(beta):
        azim = -60
    _plot_surface(ax, i1, beta, psiq,
                  (u'I1/A', u'Beta/°', u'Psi q/Vs'),
                  azim=azim)


def idq_torque(id, iq, torque, ax=0):
    """creates a surface plot of torque vs id, iq"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    unit = 'Nm'
    scale = 1
    if np.min(torque) < -9.9e3 or np.max(torque) > 9.9e3:
        scale = 1e-3
        unit = 'kNm'
    _plot_surface(ax, id, iq, scale*np.asarray(torque),
                  (u'Id/A', u'Iq/A', u'Torque/{}'.format(unit)),
                  azim=-60)
    return ax


def idq_psid(id, iq, psid, ax=0):
    """creates a surface plot of psid vs id, iq"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, id, iq, psid,
                  (u'Id/A', u'Iq/A', u'Psi d/Vs'),
                  azim=210)


def idq_psiq(id, iq, psiq, ax=0):
    """creates a surface plot of psiq vs id, iq"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, id, iq, psiq,
                  (u'Id/A', u'Iq/A', u'Psi q/Vs'),
                  azim=210)


def idq_psim(id, iq, psim, ax=0):
    """creates a surface plot of psim vs. id, iq"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, id, iq, psim,
                  (u'Id/A', u'Iq/A', u'Psi m [Vs]'),
                  azim=120)


def idq_ld(id, iq, ld, ax=0):
    """creates a surface plot of ld vs. id, iq"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, id, iq, np.asarray(ld)*1e3,
                  (u'Id/A', u'Iq/A', u'L d/mH'),
                  azim=120)


def idq_lq(id, iq, lq, ax=0):
    """creates a surface plot of lq vs. id, iq"""
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    _plot_surface(ax, id, iq, np.asarray(lq)*1e3,
                  (u'Id/A', u'Iq/A', u'L q/mH'),
                  azim=120)


def ldlq(bch):
    """creates the surface plots of a BCH reader object
    with a ld-lq identification"""
    beta = bch.ldq['beta']
    i1 = bch.ldq['i1']
    torque = bch.ldq['torque']
    ld = np.array(bch.ldq['ld'])
    lq = np.array(bch.ldq['lq'])
    psid = bch.ldq['psid']
    psiq = bch.ldq['psiq']

    rows = 3
    fig = plt.figure(figsize=(10, 4*rows))
    fig.suptitle('Ld-Lq Identification {}'.format(bch.filename), fontsize=16)
    fig.add_subplot(rows, 2, 1, projection='3d')
    i1beta_torque(i1, beta, torque)

    fig.add_subplot(rows, 2, 2, projection='3d')
    i1beta_psid(i1, beta, psid)

    fig.add_subplot(rows, 2, 3, projection='3d')
    i1beta_psiq(i1, beta, psiq)

    fig.add_subplot(rows, 2, 4, projection='3d')
    try:
        i1beta_psim(i1, beta, bch.ldq['psim'])
    except:
        i1beta_up(i1, beta, bch.ldq['up'])

    fig.add_subplot(rows, 2, 5, projection='3d')
    i1beta_ld(i1, beta, ld)

    fig.add_subplot(rows, 2, 6, projection='3d')
    i1beta_lq(i1, beta, lq)


def psidq(bch):
    """creates the surface plots of a BCH reader object
    with a psid-psiq identification"""
    id = bch.psidq['id']
    iq = bch.psidq['iq']
    torque = bch.psidq['torque']
    ld = np.array(bch.psidq_ldq['ld'])
    lq = np.array(bch.psidq_ldq['lq'])
    psim = bch.psidq_ldq['psim']
    psid = bch.psidq['psid']
    psiq = bch.psidq['psiq']

    rows = 3
    fig = plt.figure(figsize=(10, 4*rows))
    fig.suptitle('Psid-Psiq Identification {}'.format(
        bch.filename), fontsize=16)

    fig.add_subplot(rows, 2, 1, projection='3d')
    idq_torque(id, iq, torque)

    fig.add_subplot(rows, 2, 2, projection='3d')
    idq_psid(id, iq, psid)

    fig.add_subplot(rows, 2, 3, projection='3d')
    idq_psiq(id, iq, psiq)

    fig.add_subplot(rows, 2, 4, projection='3d')
    idq_psim(id, iq, psim)

    fig.add_subplot(rows, 2, 5, projection='3d')
    idq_ld(id, iq, ld)

    fig.add_subplot(rows, 2, 6, projection='3d')
    idq_lq(id, iq, lq)



def main():
    import io
    import sys
    import argparse
    from ..__init__ import __version__
    from femagtools.bch import Reader

    argparser = argparse.ArgumentParser(
        description='Read BCH/BATCH/PLT file and create a plot')
    argparser.add_argument('filename',
                           help='name of BCH/BATCH/PLT file')
    argparser.add_argument(
        "--version",
        "-v",
        action="version",
        version="%(prog)s {}, Python {}".format(__version__, sys.version),
        help="display version information",
    )
    args = argparser.parse_args()
    if not matplotlibversion:
        sys.exit(0)
    if not args.filename:
        sys.exit(0)

    ext = args.filename.split('.')[-1].upper()
    if ext.startswith('MC'):
        import femagtools.mcv
        from femagtools.plot.mcv import mcv_hbj, mcv_muer
        mcv = femagtools.mcv.read(sys.argv[1])

        if mcv['mc1_type'] in (femagtools.mcv.MAGCRV, femagtools.mcv.ORIENT_CRV):
            ncols = 2
        else:  # Permanent Magnet
            ncols = 1

        fig, ax = plt.subplots(nrows=1, ncols=ncols, figsize=(10, 6))
        if ncols > 1:
            plt.subplot(1, 2, 1)
            mcv_hbj(mcv)
            plt.subplot(1, 2, 2)
            mcv_muer(mcv)
        else:
            mcv_hbj(mcv, log=False)

        fig.tight_layout()
        fig.subplots_adjust(top=0.94)
        plt.show()
        return

    if ext.startswith('PLT'):
        import femagtools.forcedens
        from femagtools.plot.forcedens import forcedens, forcedens_fft
        fdens = femagtools.forcedens.read(args.filename)
        cols = 1
        rows = 2
        fig, ax = plt.subplots(nrows=rows, ncols=cols,
                               figsize=(10, 10*rows))
        title = '{}, Rotor position {}'.format(
            fdens.title, fdens.positions[0]['position'])
        pos = fdens.positions[0]['X']
        FT_FN = (fdens.positions[0]['FT'],
                 fdens.positions[0]['FN'])
        plt.subplot(rows, cols, 1)
        forcedens(title, pos, FT_FN)

        title = 'Force Density Harmonics'
        plt.subplot(rows, cols, 2)
        forcedens_fft(title, fdens)

        # fig.tight_layout(h_pad=3.5)
        # if title:
        #    fig.subplots_adjust(top=0.92)
        plt.show()
        return

    bchresults = Reader()
    with io.open(args.filename, encoding='latin1', errors='ignore') as f:
        bchresults.read(f.readlines())

    if (bchresults.type.lower().find(
        'pm-synchronous-motor simulation') >= 0 or
        bchresults.type.lower().find(
            'permanet-magnet-synchronous-motor') >= 0 or
        bchresults.type.lower().find(
            'simulation pm/universal-motor') >= 0):
        pmrelsim(bchresults, bchresults.filename)
    elif bchresults.type.lower().find(
            'multiple calculation of forces and flux') >= 0:
        multcal(bchresults, bchresults.filename)
    elif bchresults.type.lower().find('cogging calculation') >= 0:
        cogging(bchresults, bchresults.filename)
    elif bchresults.type.lower().find('ld-lq-identification') >= 0:
        ldlq(bchresults)
    elif bchresults.type.lower().find('psid-psiq-identification') >= 0:
        psidq(bchresults)
    elif bchresults.type.lower().find('fast_torque calculation') >= 0:
        fasttorque(bchresults)
    elif bchresults.type.lower().find('transient sc') >= 0:
        transientsc(bchresults, bchresults.filename)
    else:
        raise ValueError("BCH type {} not yet supported".format(
            bchresults.type))
    plt.show()


def eigenmode(reigen, num_modes=12):
    """plot eigenmode"""
    cols = 4
    rows = int((num_modes - num_modes%cols)/cols + 1)
    fig, ax = plt.subplots(rows, cols)
    ctr = 0
    for i in range(rows):
        for j in range(cols):
            if ctr < num_modes:
                ax[i, j].imshow(reigen[0][ctr])
                ax[i, j].set_title(f'Freq: {reigen[1][ctr]:.2f} Hz')
                ax[i, j].axis('off')
            else:
                ax[i, j].remove()
                ax[i, j] = None
            ctr+=1
    plt.tight_layout()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    main()
