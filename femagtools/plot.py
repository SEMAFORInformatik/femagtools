"""
    femagtools.plot
    ~~~~~~~~~~~~~~~

    Creating plots

    :copyright: 2017 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import os
import matplotlib.pyplot as pl
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as ip


def plot_surface(ax, x, y, z, labels, azim=None):
    """xticks, yticks, zticks):"""
    ax.tick_params(axis='both', which='major', pad=-3)
    if azim is not None:
        ax.azim = azim
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y, z,
                    rstride=1, cstride=1, cmap=cm.jet,
                    linewidth=0, antialiased=True,
                    edgecolor=(0, 0, 0, 0))

    #ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    #ax.set_zticks(zticks)

    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(labels[2])
    
    #pl.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)


def phasor_plot(up, i1, beta, r1, xd, xq, file=None):
    """creates a phasor plot
    up: internal voltage
    i1: current
    beta: angle i1 vs up [deg]
    r1: resistance
    xd: reactance in direct axis
    xq: reactance in quadrature axis"""

    i1d, i1q = (i1*np.sin(beta/180*np.pi), i1*np.cos(beta/180*np.pi))
    uxd, uxq = ((r1*i1d - xq*i1q)/up, (r1*i1q + xd*i1d)/up)

    u1d, u1q = (uxd, 1+uxq)
    i1d, i1q = (i1d/i1, i1q/i1)

    hw = 0.05
    ax = pl.axes(aspect='equal')
    ax.arrow(0, 0, 0, 1, color='k', head_width=hw,
             length_includes_head=True)
    ax.text(0.08, 0.9, r'$U_p$', fontsize=18)
    ax.arrow(0, 0, u1d, u1q, color='r', head_width=hw,
             length_includes_head=True)
    ax.text(u1d, 0.75*u1q, r'$U_1$', fontsize=18)
    ax.arrow(0, 1, uxd, 0, color='g', head_width=hw,
             length_includes_head=True)
    ax.arrow(uxd, 1, 0, uxq, color='g', head_width=hw,
             length_includes_head=True)
    ax.arrow(0, 0, i1d, i1q, color='b', head_width=hw,
             length_includes_head=True)
    ax.text(1.15*i1d, 0.72*i1q, r'$I_1$', fontsize=18)

    xmin, xmax = (min(0, uxd, i1d), max(0, i1d, uxd))
    ymin, ymax = (min(0, i1q, uxq), max(1, i1q))
    pl.axis([xmin-0.1, xmax+0.1, ymin-0.1, ymax+0.1])
    pl.grid(True)
    if file:
        pl.savefig(file)
    else:
        pl.show()
        

def pmrelsim_plot(bch, title='', file=None):
    """creates a plot of a PM/Rel motor simulation"""
    torque = bch.torque[-1]
    torque_fft = bch.torque_fft[-1]
    flux = [bch.flux[k][-1] for k in bch.flux]
    flux_fft = bch.flux_fft['1'][-1]
    k = 20

    alpha = np.linspace(0, max(torque['angle']),
                        k*len(torque['torque']))
    f = ip.interp1d(torque['angle'], torque['torque'], kind='cubic')
    torque['angle_fit'] = alpha
    torque['torque_spline'] = f(alpha)

    cases = [
        dict(
            x=(torque['angle'],
               torque['angle_fit']),
            y=(torque['torque'],
               torque['torque_spline']),
            title='Torque [Nm]'),
        dict(
            x=torque_fft['order'],
            y=torque_fft['torque'],
            title='Torque Harmonics [Nm]'),

        dict(
            x=torque['angle'],
            y=torque['force_x'],
            title='Force Fx [N]'),
        dict(
            x=torque['angle'],
            y=torque['force_y'],
            title='Force Fy [N]'),

        dict(
            x=(flux[0]['displ'],
               flux[1]['displ'],
               flux[2]['displ']),
            y=(flux[0]['flux_k'],
               flux[1]['flux_k'],
               flux[2]['flux_k']),
            title='Winding Flux [Vs]'),
        dict(
            x=(flux[0]['displ'],
               flux[1]['displ'],
               flux[2]['displ']),
            y=(flux[0]['current_k'],
               flux[1]['current_k'],
               flux[2]['current_k']),
            title='Winding Current [A]'),

        dict(
            x=flux[0]['displ'],
            y=flux[0]['voltage_dpsi'],
            title='Internal Voltage [V]'),
        dict(
            x=flux_fft['order'],
            y=flux_fft['voltage'],
            title='Internal Voltage Harmonics [V]')]

    if len(bch.flux['1']) > 1:
        cases += [
            dict(
                x=bch.flux['1'][0]['displ'],
                y=bch.flux['1'][0]['voltage_dpsi'],
                title='No Load Voltage [V]'),
            dict(
                x=bch.flux_fft['1'][0]['order'],
                y=bch.flux_fft['1'][0]['voltage'],
                title='No Load Voltage Harmonics [V]')]

    cols = 2
    rows = len(cases) // cols
    htitle = 1.5 if title else 0
    fig, ax = pl.subplots(nrows=rows, ncols=cols,
                          figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)
    
    ax[0, 0].set_title(cases[0]['title'])
    ax[0, 0].grid(True)

    ax[0, 0].plot(cases[0]['x'][0], cases[0]['y'][0], 'go')
    ax[0, 0].plot(cases[0]['x'][1], cases[0]['y'][1])
    ax[0, 0].set_ylim(bottom=0)

    ax[0, 1].set_title(cases[1]['title'])
    ax[0, 1].grid(True)

    bw = 2.5E-2*max(cases[1]['x'])
    ax[0, 1].bar(cases[1]['x'], cases[1]['y'], width=bw, align='center')
    ax[0, 1].set_xlim(left=0)

    ax[1, 0].set_title(cases[2]['title'])
    ax[1, 0].grid(True)

    ax[1, 0].plot(cases[2]['x'], cases[2]['y'])
    ax[1, 0].set_ylim(bottom=0)

    ax[1, 1].set_title(cases[3]['title'])
    ax[1, 1].grid(True)

    ax[1, 1].plot(cases[3]['x'], cases[3]['y'])
    ax[1, 1].set_ylim(bottom=0)

    ax[2, 0].set_title(cases[4]['title'])
    ax[2, 0].grid(True)

    ax[2, 0].plot(cases[4]['x'][0], cases[4]['y'][0])
    ax[2, 0].plot(cases[4]['x'][1], cases[4]['y'][1])
    ax[2, 0].plot(cases[4]['x'][2], cases[4]['y'][2])

    ax[2, 1].set_title(cases[5]['title'])
    ax[2, 1].grid(True)

    ax[2, 1].plot(cases[5]['x'][0], cases[5]['y'][0])
    ax[2, 1].plot(cases[5]['x'][1], cases[5]['y'][1])
    ax[2, 1].plot(cases[5]['x'][2], cases[5]['y'][2])

    ax[3, 0].set_title('Internal Voltage [V]')
    ax[3, 0].grid(True)

    ax[3, 0].plot(cases[6]['x'], cases[6]['y'])

    ax[3, 1].set_title(cases[7]['title'])
    ax[3, 1].grid(True)
    bw = 2.5E-2*max(cases[7]['x'])
    ax[3, 1].bar(cases[7]['x'], cases[7]['y'], width=bw, align='center')

    if len(cases) > 8:
        ax[4, 0].set_title('No Load Voltage [V]')
        ax[4, 0].grid(True)

        ax[4, 0].plot(cases[8]['x'], cases[8]['y'])

        ax[4, 1].set_title(cases[9]['title'])
        ax[4, 1].grid(True)

        bw = 2.5E-2*max(cases[9]['x'])
        ax[4, 1].bar(cases[9]['x'], cases[9]['y'], width=bw, align='center')

    if title:
        fig.tight_layout(h_pad=2.5)
        fig.subplots_adjust(top=0.9)
    else:
        fig.tight_layout(h_pad=2.5)
       
    if file:
        pl.savefig(file)
    else:
        pl.show()
    pl.close()


def cogging_plot(bch, file=None):
    torque = bch.torque[0]
    torque_fft = bch.torque_fft[0]
    flux = [bch.flux[k][0] for k in bch.flux]
    flux_fft = bch.flux_fft['1'][0]
    k = 20

    alpha = np.linspace(0, max(torque['angle']),
                        k*len(torque['torque']))
    f = ip.interp1d(torque['angle'], torque['torque'], kind='cubic')
    torque['angle_fit'] = alpha
    torque['torque_spline'] = f(alpha)

    cases = [
        dict(
            x=(torque['angle'],
               torque['angle_fit']),
            y=(torque['torque'],
               torque['torque_spline']),
            title='Torque [Nm]'),
        dict(
            x=torque_fft['order'],
            y=torque_fft['torque'],
            title='Torque Harmonics [Nm]'),

        dict(
            x=torque['angle'],
            y=torque['force_x'],
            title='Force Fx [N]'),
        dict(
            x=torque['angle'],
            y=torque['force_y'],
            title='Force Fy [N]'),
        
        dict(
            x=flux[0]['displ'],
            y=(flux[0]['voltage_dpsi'],
               flux[1]['voltage_dpsi'],
               flux[2]['voltage_dpsi']),
            title='Voltage [V]'),
        dict(
            x=flux_fft['order'],
            y=flux_fft['voltage'],
            title='Voltage Harmonics [V]')]

    cols = 2
    rows = len(cases) // cols
    fig1, ax = pl.subplots(ncols=cols, nrows=rows, figsize=(10, 3*rows))

    ax[0, 0].set_title(cases[0]['title'])
    ax[0, 0].grid(True)

    ax[0, 0].plot(cases[0]['x'][0], cases[0]['y'][0], 'go')
    ax[0, 0].plot(cases[0]['x'][1], cases[0]['y'][1])

    ax[0, 1].set_title(cases[1]['title'])
    ax[0, 1].grid(True)

    bw = 2.5E-2*max(cases[1]['x'])
    ax[0, 1].bar(cases[1]['x'], cases[1]['y'], width=bw, align='center')
    ax[0, 1].set_xlim(left=0)
    
    ax[1, 0].set_title(cases[0]['title'])
    ax[1, 0].grid(True)

    ax[1, 0].set_title(cases[2]['title'])
    ax[1, 0].grid(True)

    ax[1, 0].plot(cases[2]['x'], cases[2]['y'])
    ax[1, 0].set_ylim(bottom=0)

    ax[1, 1].set_title(cases[3]['title'])
    ax[1, 1].grid(True)

    ax[1, 1].plot(cases[3]['x'], cases[3]['y'])
    ax[1, 1].set_ylim(bottom=0)

    ax[2, 0].set_title(cases[4]['title'])
    ax[2, 0].grid(True)

    ax[2, 0].plot(cases[4]['x'], cases[4]['y'][0])
    ax[2, 0].plot(cases[4]['x'], cases[4]['y'][1])
    ax[2, 0].plot(cases[4]['x'], cases[4]['y'][2])

    ax[2, 1].set_title(cases[5]['title'])
    ax[2, 1].grid(True)

    bw = 2.5E-2*max(cases[5]['x'])
    ax[2, 1].bar(cases[5]['x'], cases[5]['y'], width=bw, align='center')
    ax[2, 1].set_xlim(left=0)

    fig1.tight_layout()
    if file:
        pl.savefig(file)
    else:
        pl.show()
    pl.close()


def ldlq_plot(bch):
    beta = bch.ldq['beta']
    i1 = bch.ldq['i1']
    torque = bch.ldq['torque']
    ld = np.array(bch.ldq['ld'])*1e3
    lq = np.array(bch.ldq['lq'])*1e3
    psim = bch.ldq['psim']
    project = os.path.basename(bch.project).split('.')[0]

    rows = 2
    fig = pl.figure(figsize=(10, 3.5*rows))
    fig.suptitle('Ld-Lq Identification {}'.format(project), fontsize=16)
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    plot_surface(ax, i1, beta, torque,
                 (u'I1 / A', u'Beta / °', u'Torque / Nm'),
                 azim=210)
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    plot_surface(ax, i1, beta, ld,
                 (u'I1 / A', u'Beta / °', u'Ld / mH'))
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    plot_surface(ax, i1, beta, lq,
                 (u'I1 / A', u'Beta / °', u'Lq / mH'))
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    plot_surface(ax, i1, beta, psim,
                 (u'I1 / A', u'Beta / °', u'Psi m /Vs'))
    pl.show()
        

def psidq_plot(bch):
    id = bch.psidq['id']
    iq = bch.psidq['iq']
    torque = bch.psidq['torque']
    #ld = np.array(bch.psidq['ld'])*1e3
    #lq = np.array(bch.psidq['lq'])*1e3
    psid = bch.psidq['psid']
    psiq = bch.psidq['psiq']
    project = bch.project

    rows = 2
    fig = pl.figure(figsize=(10, 3.5*rows))
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    plot_surface(ax, id, iq, torque,
                 (u'Id /A', u'Iq /A', u'Torque /Nm'),
                 azim=210)
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    plot_surface(ax, id, iq, psid,
                 (u'Id /A', u'Iq /A', u'Psi d /Vs'))
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    plot_surface(ax, id, iq, psiq,
                 (u'Id /A', u'Iq /A', u'Psi q /Vs'))
    #ax = fig.add_subplot(2, 2, 4, projection='3d')
    #plot_surface(ax, id, iq, psim,
    #             (u'Id /A', u'Iq /A', u'Psi m [Vs]'))
    pl.show()

if __name__ == "__main__":
    import io
    import sys
    from femagtools.bch import Reader
    for filename in sys.argv[1:]:
        bchresults = Reader()
        with io.open(filename, encoding='latin1', errors='ignore') as f:
            bchresults.read(f.readlines())

        if bchresults.type == 'Fast PM-Synchronous-Motor Simulation':
            pmrelsim_plot(bchresults)
        elif bchresults.type == 'Fast cogging calculation OF FORCES AND FLUX':
            cogging_plot(bchresults)
        elif bchresults.type == 'Fast LD-LQ-Identification':
            ldlq_plot(bchresults)
        elif bchresults.type == 'Fast Psid-Psiq-Identification':
            psidq_plot(bchresults)
    
    
