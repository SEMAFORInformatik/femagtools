# -*- coding: utf-8 -*-
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
    u1 = np.sqrt(u1d**2 + u1q**2)*up
    i1d, i1q = (i1d/i1, i1q/i1)

    hw = 0.05
    ax = pl.gca()
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    
    ax.set_title(
        r'$U_1$={0} V, $I_1$={1} A, $U_p$={2} V'.format(
            round(u1, 1), round(i1, 1), round(up, 1)), fontsize=14)
    ax.set_aspect('equal', adjustable='box')
    ax.arrow(0, 0, 0, 1, color='k', head_width=hw,
             length_includes_head=True)
    ax.text(0.08, 0.9, r'$U_p$', fontsize=18)
    ax.arrow(0, 0, u1d, u1q, color='r', head_width=hw,
             length_includes_head=True)
    ax.text(u1d, 0.75*u1q, r'$U_1$', fontsize=18)
    ax.arrow(0, 1, uxd, 0, color='g', head_width=hw,
             length_includes_head=True)
    if abs(uxq) > 1e-3:
        ax.arrow(uxd, 1, 0, uxq, color='g', head_width=hw,
                 length_includes_head=True)
    ax.arrow(0, 0, i1d, i1q, color='b', head_width=hw,
             length_includes_head=True)
    ax.text(1.15*i1d, 0.72*i1q, r'$I_1$', fontsize=18)

    xmin, xmax = (min(0, uxd, i1d), max(0, i1d, uxd))
    ymin, ymax = (min(0, i1q, uxq), max(1, i1q))
    ax.set_xlim([xmin-0.1, xmax+0.1])
    ax.set_ylim([ymin-0.1, ymax+0.1])
    ax.grid(True)
        
def torque_plot(torque):
    """plot torque vs position"""
    k = 20
    alpha = np.linspace(0, max(torque['angle']),
                        k*len(torque['torque']))
    f = ip.interp1d(torque['angle'], torque['torque'], kind='cubic')
    ax = pl.gca()
    ax.set_title('Torque / Nm')
    ax.grid(True)
    ax.plot(torque['angle'], torque['torque'], 'go')
    ax.plot(alpha, f(alpha))
    if min(torque['torque']) > 0:
        ax.set_ylim(bottom=0)
    
def torque_fft_plot(torque_fft):
    """plot torque harmonics"""
    ax = pl.gca()
    ax.set_title('Torque Harmonics [Nm]')
    ax.grid(True)
    bw = 2.5E-2*max(torque_fft['order'])
    ax.bar(torque_fft['order'], torque_fft['torque'], width=bw, align='center')
    ax.set_xlim(left=-bw/2)

def force_x_plot(torque):
    """plot force x vs position"""
    x = torque['angle']
    y = torque['force_x']
    ax = pl.gca()
    ax.set_title('Force Fx [N]')
    ax.grid(True)
    ax.plot(x, y)
    ax.set_ylim(bottom=0)

def force_y_plot(torque):
    """plot force y vs position"""
    x = torque['angle']
    y = torque['force_y']
    ax = pl.gca()
    ax.set_title('Force Fy [N]')
    ax.grid(True)
    ax.plot(x, y)
    ax.set_ylim(bottom=0)

def winding_flux_plot(flux):
    ax = pl.gca()
    ax.set_title('Winding Flux [Vs]')
    ax.grid(True)
    ax.plot(flux[0]['displ'], flux[0]['flux_k'])
    ax.plot(flux[1]['displ'], flux[1]['flux_k'])
    ax.plot(flux[2]['displ'], flux[2]['flux_k'])

def winding_current_plot(flux):
    ax = pl.gca()
    ax.set_title('Winding Currents [A]')
    ax.grid(True)
    ax.plot(flux[0]['displ'], flux[0]['current_k'])
    ax.plot(flux[1]['displ'], flux[1]['current_k'])
    ax.plot(flux[2]['displ'], flux[2]['current_k'])

def voltage_plot(title, pos, voltage):
    ax = pl.gca()
    ax.set_title(title)
    ax.grid(True)
    ax.plot(pos, voltage)
    
def voltage_fft_plot(title, order, voltage):
    ax = pl.gca()
    ax.set_title('{} [V]'.format(title))
    ax.grid(True)
    bw = 2.5E-2*max(order)
    ax.bar(order, voltage, width=bw, align='center')


def pmrelsim_plot(bch, title='', file=None):
    """creates a plot of a PM/Rel motor simulation"""
    cols = 2
    rows = 4
    if len(bch.flux['1']) > 1:
        rows += 1
    htitle = 1.5 if title else 0
    fig, ax = pl.subplots(nrows=rows, ncols=cols,
                          figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)
    
    pl.subplot(rows, cols, 1)
    torque_plot(bch.torque[-1])
    pl.subplot(rows, cols, 2)
    torque_fft_plot(bch.torque_fft[-1])
    pl.subplot(rows, cols, 3)
    force_x_plot(bch.torque[-1])
    pl.subplot(rows, cols, 4)
    force_y_plot(bch.torque[-1])
    pl.subplot(rows, cols, 5)
    flux = [bch.flux[k][-1] for k in bch.flux]
    winding_flux_plot(flux)
    pl.subplot(rows, cols, 6)
    winding_current_plot(flux)
    pl.subplot(rows, cols, 7)
    voltage_plot('Internal Voltage [V]',
                 bch.flux['1'][-1]['displ'],
                 bch.flux['1'][-1]['voltage_dpsi'])
    pl.subplot(rows, cols, 8)
    voltage_fft_plot('Internal Voltage Harmonics [V]',
                     bch.flux_fft['1'][-1]['order'],
                     bch.flux_fft['1'][-1]['voltage'])
                              
    if len(bch.flux['1']) > 1:
        pl.subplot(rows, cols, 9)
        voltage_plot('No Load Voltage [V]',
                     bch.flux['1'][0]['displ'],
                     bch.flux['1'][0]['voltage_dpsi'])
        pl.subplot(rows, cols, 10)
        voltage_fft_plot('No Load Voltage Harmonics [V]',
                         bch.flux_fft['1'][-1]['order'],
                         bch.flux_fft['1'][-1]['voltage'])

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


def cogging_plot(bch, title='', file=None):
    """creates a cogging plot"""
    flux = [bch.flux[k][0] for k in bch.flux]

    cols = 2
    rows = 4
    htitle = 1.5 if title else 0
    fig, ax = pl.subplots(nrows=rows, ncols=cols,
                          figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)
    
    pl.subplot(rows, cols, 1)
    torque_plot(bch.torque[0])
    pl.subplot(rows, cols, 2)
    torque_fft_plot(bch.torque_fft[0])
    pl.subplot(rows, cols, 3)
    force_x_plot(bch.torque[0])
    pl.subplot(rows, cols, 4)
    force_y_plot(bch.torque[0])
    pl.subplot(rows, cols, 5)
    flux = [bch.flux[k][-1] for k in bch.flux]
    winding_flux_plot(flux)
    pl.subplot(rows, cols, 6)
    winding_current_plot(flux)
    pl.subplot(rows, cols, 7)
    voltage_plot('Voltage [V]',
                 bch.flux['1'][0]['displ'],
                 bch.flux['1'][0]['voltage_dpsi'])
    pl.subplot(rows, cols, 8)
    voltage_fft_plot('Voltage Harmonics [V]',
                     bch.flux_fft['1'][0]['order'],
                     bch.flux_fft['1'][0]['voltage'])

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
    
    
