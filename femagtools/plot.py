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


def _create_3d_axis():
    """creates a subplot with 3d projection if one does not already exist"""
    from matplotlib.projections import get_projection_class
    from matplotlib import _pylab_helpers

    create_axis = True
    if _pylab_helpers.Gcf.get_active() is not None:
        if isinstance(pl.gca(), get_projection_class('3d')):
            create_axis = False
    if create_axis:
        pl.figure()
        pl.subplot(111, projection='3d')

        
def _plot_surface(ax, x, y, z, labels, azim=None):
    """helper function for surface plots"""
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


def __phasor_plot(up, idq, uxdq):
    uxd = uxdq[0]/up
    uxq = uxdq[1]/up
    u1d, u1q = (uxd, 1+uxq)
    u1 = np.sqrt(u1d**2 + u1q**2)*up
    i1 = np.linalg.norm(idq)
    i1d, i1q = (idq[0]/i1, idq[1]/i1)

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


def i1beta_phasor_plot(up, i1, beta, r1, xd, xq):
    """creates a phasor plot
    up: internal voltage
    i1: current
    beta: angle i1 vs up [deg]
    r1: resistance
    xd: reactance in direct axis
    xq: reactance in quadrature axis"""

    i1d, i1q = (i1*np.sin(beta/180*np.pi), i1*np.cos(beta/180*np.pi))
    uxdq = ((r1*i1d - xq*i1q), (r1*i1q + xd*i1d))
    __phasor_plot(up, (i1d, i1q), uxdq)


def iqd_phasor_plot(up, iqd, uqd):
    """creates a phasor plot
    up: internal voltage
    iqd: current
    uqd: terminal voltage"""

    uxdq = (uqd[1]/np.sqrt(2), (uqd[0]/np.sqrt(2)-up))
    __phasor_plot(up, (iqd[1]/np.sqrt(2), iqd[0]/np.sqrt(2)), uxdq)


def phasor_plot(bch):
    """create phasor plot from bch"""
    f1 = bch.machine['p']*bch.dqPar['speed']
    w1 = 2*np.pi*f1
    xd = w1*bch.dqPar['ld'][1]
    xq = w1*bch.dqPar['lq'][1]
    r1 = bch.machine['r1']
    i1beta_phasor_plot(bch.dqPar['up0'],
                       bch.dqPar['i1'][-1], bch.dqPar['beta'][1],
                       r1, xd, xq)
    

def airgap_plot(airgap):
    """creates plot of flux density in airgap"""
    pl.title('Airgap Induction [T]')
    pl.plot(airgap['pos'], airgap['B'])
    pl.plot(airgap['pos'], airgap['B_fft'])
    pl.xlabel('Position/°')
    pl.grid()

    
def torque_plot(pos, torque):
    """creates plot from torque vs position"""
    k = 20
    alpha = np.linspace(0, pos[-1],
                        k*len(torque))
    f = ip.interp1d(pos, torque, kind='cubic')
    unit = 'Nm'
    scale = 1
    if min(torque) < -9.9e3 or max(torque) > 9.9e3:
        scale = 1e-3
        unit = 'kNm'
    ax = pl.gca()
    ax.set_title('Torque / {}'.format(unit))
    ax.grid(True)
    ax.plot(pos, [scale*t for t in torque], 'go')
    ax.plot(alpha, scale*f(alpha))
    if min(torque) > 0 and max(torque) > 0:
        ax.set_ylim(bottom=0)
    elif min(torque) < 0 and max(torque) < 0:
        ax.set_ylim(top=0)
    

def torque_fft_plot(order, torque):
    """plot torque harmonics"""
    unit = 'Nm'
    scale = 1
    if min(torque) < -9.9e3 or max(torque) > 9.9e3:
        scale = 1e-3
        unit = 'kNm'
    ax = pl.gca()
    ax.set_title('Torque Harmonics / {}'.format(unit))
    ax.grid(True)
    bw = 2.5E-2*max(order)
    ax.bar(order, [scale*t for t in torque], width=bw, align='center')
    ax.set_xlim(left=-bw/2)


def force_plot(title, pos, force):
    """plot force vs position"""
    unit = 'N'
    scale = 1
    if min(force) < -9.9e3 or max(force) > 9.9e3:
        scale = 1e-3
        unit = 'kN'
    ax = pl.gca()
    ax.set_title('{} / {}'.format(title, unit))
    ax.grid(True)
    ax.plot(pos, [scale*f for f in force])
    ax.set_ylim(bottom=0)


def winding_flux_plot(pos, flux):
    ax = pl.gca()
    ax.set_title('Winding Flux / Vs')
    ax.grid(True)
    ax.plot(pos[0], flux[0])
    ax.plot(pos[1], flux[1])
    ax.plot(pos[2], flux[2])


def winding_current_plot(pos, current):
    ax = pl.gca()
    ax.set_title('Winding Currents / A')
    ax.grid(True)
    ax.plot(pos[0], current[0])
    ax.plot(pos[1], current[1])
    ax.plot(pos[2], current[2])


def voltage_plot(title, pos, voltage):
    ax = pl.gca()
    ax.set_title('{} / V'.format(title))
    ax.grid(True)
    ax.plot(pos, voltage)
    

def voltage_fft_plot(title, order, voltage):
    ax = pl.gca()
    ax.set_title('{} / V'.format(title))
    ax.grid(True)
    bw = 2.5E-2*max(order)
    ax.bar(order, voltage, width=bw, align='center')


def pmrelsim_plot(bch, title=''):
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
    torque_plot(bch.torque[-1]['angle'], bch.torque[-1]['torque'])
    pl.subplot(rows, cols, 2)
    torque_fft_plot(bch.torque_fft[-1]['order'], bch.torque_fft[-1]['torque'])
    pl.subplot(rows, cols, 3)
    force_plot('Force Fx',
               bch.torque[-1]['angle'], bch.torque[-1]['force_x'])
    pl.subplot(rows, cols, 4)
    force_plot('Force Fy',
               bch.torque[-1]['angle'], bch.torque[-1]['force_y'])
    pl.subplot(rows, cols, 5)
    flux = [bch.flux[k][-1] for k in bch.flux]
    pos = [f['displ'] for f in flux]
    winding_flux_plot(pos,
                      (flux[0]['flux_k'],
                       flux[1]['flux_k'],
                       flux[2]['flux_k']))
    pl.subplot(rows, cols, 6)
    winding_current_plot(pos,
                         (flux[0]['current_k'],
                          flux[1]['current_k'],
                          flux[2]['current_k']))
    pl.subplot(rows, cols, 7)
    voltage_plot('Internal Voltage',
                 bch.flux['1'][-1]['displ'],
                 bch.flux['1'][-1]['voltage_dpsi'])
    pl.subplot(rows, cols, 8)
    voltage_fft_plot('Internal Voltage Harmonics',
                     bch.flux_fft['1'][-1]['order'],
                     bch.flux_fft['1'][-1]['voltage'])
                              
    if len(bch.flux['1']) > 1:
        pl.subplot(rows, cols, 9)
        voltage_plot('No Load Voltage',
                     bch.flux['1'][0]['displ'],
                     bch.flux['1'][0]['voltage_dpsi'])
        pl.subplot(rows, cols, 10)
        voltage_fft_plot('No Load Voltage Harmonics',
                         bch.flux_fft['1'][-1]['order'],
                         bch.flux_fft['1'][-1]['voltage'])

    fig.tight_layout(h_pad=3.5)
    if title:
        fig.subplots_adjust(top=0.92)
    

def cogging_plot(bch, title=''):
    """creates a cogging plot"""
    cols = 2
    rows = 3
    htitle = 1.5 if title else 0
    fig, ax = pl.subplots(nrows=rows, ncols=cols,
                          figsize=(10, 3*rows + htitle))
    if title:
        fig.suptitle(title, fontsize=16)
    
    pl.subplot(rows, cols, 1)
    torque_plot(bch.torque[0]['angle'], bch.torque[0]['torque'])
    pl.subplot(rows, cols, 2)
    torque_fft_plot(bch.torque_fft[0]['order'], bch.torque_fft[0]['torque'])
    pl.subplot(rows, cols, 3)
    force_plot('Force Fx',
               bch.torque[0]['angle'], bch.torque[0]['force_x'])
    pl.subplot(rows, cols, 4)
    force_plot('Force Fy',
               bch.torque[0]['angle'], bch.torque[0]['force_y'])
    pl.subplot(rows, cols, 5)
    voltage_plot('Voltage',
                 bch.flux['1'][0]['displ'],
                 bch.flux['1'][0]['voltage_dpsi'])
    pl.subplot(rows, cols, 6)
    voltage_fft_plot('Voltage Harmonics',
                     bch.flux_fft['1'][0]['order'],
                     bch.flux_fft['1'][0]['voltage'])

    fig.tight_layout(h_pad=2)
    if title:
        fig.subplots_adjust(top=0.92)
       

def i1beta_torque_plot(i1, beta, torque):
    """creates a surface plot of torque vs i1, beta"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, i1, beta, torque,
                  (u'I1/A', u'Beta/°', u'Torque/Nm'),
                  azim=210)


def i1beta_ld_plot(i1, beta, ld):
    """creates a surface plot of ld vs i1, beta"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, i1, beta, ld,
                  (u'I1/A', u'Beta/°', u'Ld/mH'),
                  azim=60)
    

def i1beta_lq_plot(i1, beta, lq):
    """creates a surface plot of ld vs i1, beta"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, i1, beta, lq,
                  (u'I1/A', u'Beta/°', u'Lq/mH'),
                  azim=60)


def i1beta_psim_plot(i1, beta, psim):
    """creates a surface plot of psim vs i1, beta"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, i1, beta, psim,
                  (u'I1/A', u'Beta/°', u'Psi m/Vs'),
                  azim=60)


def i1beta_psid_plot(i1, beta, psid):
    """creates a surface plot of psid vs i1, beta"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, i1, beta, psid,
                  (u'I1/A', u'Beta/°', u'Psi d/Vs'),
                  azim=-60)


def i1beta_psiq_plot(i1, beta, psiq):
    """creates a surface plot of psiq vs i1, beta"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, i1, beta, psiq,
                  (u'I1/A', u'Beta/°', u'Psi q/Vs'),
                  azim=210)


def idq_torque_plot(id, iq, torque):
    """creates a surface plot of torque vs id, iq"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, id, iq, torque,
                  (u'Id/A', u'Iq/A', u'Torque/Nm'),
                  azim=-60)


def idq_psid_plot(id, iq, psid):
    """creates a surface plot of psid vs id, iq"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, id, iq, psid,
                  (u'Id/A', u'Iq/A', u'Psi d/Vs'),
                  azim=210)


def idq_psiq_plot(id, iq, psiq):
    """creates a surface plot of psiq vs id, iq"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, id, iq, psiq,
                  (u'Id/A', u'Iq/A', u'Psi q/Vs'),
                  azim=210)


def idq_psim_plot(id, iq, psim):
    """creates a surface plot of psim vs. id, iq"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, id, iq, psim,
                  (u'Id/A', u'Iq/A', u'Psi m [Vs]'),
                  azim=120)
    

def idq_ld_plot(id, iq, ld):
    """creates a surface plot of ld vs. id, iq"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, id, iq, ld,
                  (u'Id/A', u'Iq/A', u'L d/mH'),
                  azim=120)
    

def idq_lq_plot(id, iq, lq):
    """creates a surface plot of lq vs. id, iq"""
    _create_3d_axis()
    ax = pl.gca()
    _plot_surface(ax, id, iq, lq,
                  (u'Id/A', u'Iq/A', u'L q/mH'),
                  azim=120)
    

def ldlq_plot(bch):
    """creates the surface plots of a BCH reader object
    with a ld-lq identification"""
    beta = bch.ldq['beta']
    i1 = bch.ldq['i1']
    torque = bch.ldq['torque']
    ld = np.array(bch.ldq['ld'])*1e3
    lq = np.array(bch.ldq['lq'])*1e3
    psid = bch.ldq['psid']
    psiq = bch.ldq['psiq']
    psim = bch.ldq['psim']

    rows = 3
    fig = pl.figure(figsize=(10, 4*rows))
    fig.suptitle('Ld-Lq Identification {}'.format(bch.filename), fontsize=16)
    fig.add_subplot(rows, 2, 1, projection='3d')
    i1beta_torque_plot(i1, beta, torque)

    fig.add_subplot(rows, 2, 2, projection='3d')
    i1beta_psid_plot(i1, beta, psid)
    
    fig.add_subplot(rows, 2, 3, projection='3d')
    i1beta_psiq_plot(i1, beta, psiq)

    fig.add_subplot(rows, 2, 4, projection='3d')
    i1beta_psim_plot(i1, beta, psim)

    fig.add_subplot(rows, 2, 5, projection='3d')
    i1beta_ld_plot(i1, beta, ld)

    fig.add_subplot(rows, 2, 6, projection='3d')
    i1beta_lq_plot(i1, beta, lq)


def psidq_plot(bch):
    """creates the surface plots of a BCH reader object
    with a psid-psiq identification"""
    id = bch.psidq['id']
    iq = bch.psidq['iq']
    torque = bch.psidq['torque']
    ld = np.array(bch.psidq_ldq['ld'])*1e3
    lq = np.array(bch.psidq_ldq['lq'])*1e3
    psim = bch.psidq_ldq['psim']
    psid = bch.psidq['psid']
    psiq = bch.psidq['psiq']

    rows = 3
    fig = pl.figure(figsize=(10, 4*rows))
    fig.suptitle('Psid-Psiq Identification {}'.format(
        bch.filename), fontsize=16)

    fig.add_subplot(rows, 2, 1, projection='3d')
    idq_torque_plot(id, iq, torque)
    
    fig.add_subplot(rows, 2, 2, projection='3d')
    idq_psid_plot(id, iq, psid)
    
    fig.add_subplot(rows, 2, 3, projection='3d')
    idq_psiq_plot(id, iq, psiq)

    fig.add_subplot(rows, 2, 4, projection='3d')
    idq_psim_plot(id, iq, psim)
    
    fig.add_subplot(rows, 2, 5, projection='3d')
    idq_ld_plot(id, iq, ld)
    
    fig.add_subplot(rows, 2, 6, projection='3d')
    idq_lq_plot(id, iq, lq)

if __name__ == "__main__":
    import io
    import sys
    from femagtools.bch import Reader
    for filename in sys.argv[1:]:
        bchresults = Reader()
        with io.open(filename, encoding='latin1', errors='ignore') as f:
            bchresults.read(f.readlines())

        if bchresults.type == 'Fast PM-Synchronous-Motor Simulation':
            pmrelsim_plot(bchresults, bchresults.filename)
        elif bchresults.type == 'Fast cogging calculation OF FORCES AND FLUX':
            cogging_plot(bchresults, bchresults.filename)
        elif bchresults.type == 'Fast LD-LQ-Identification':
            ldlq_plot(bchresults)
        elif bchresults.type == 'Fast Psid-Psiq-Identification':
            psidq_plot(bchresults)
        else:
            raise ValueError("BCH type {} not yet supported".format(
                bchresults.type))
        pl.show()
