"""
    femagtools.plot.phasor
    ~~~~~~~~~~~~~~~~~~~~~~

    Creating phasor plots


"""
import numpy as np
import matplotlib.pyplot as plt


def __phasor_plot(ax, up, idq, uxdq):
    uref = max(up, uxdq[0])
    uxd = uxdq[0]/uref
    uxq = uxdq[1]/uref
    u1d, u1q = (uxd, 1+uxq)
    u1 = np.sqrt(u1d**2 + u1q**2)*uref
    i1 = np.linalg.norm(idq)
    i1d, i1q = (idq[0]/i1, idq[1]/i1)

    qhw = 6   # width arrow head
    qhl = 15  # length arrow head
    qlw = 2   # line width
    qts = 10  # textsize
    # Length of the Current adjust to Ud: Initally 0.9, Maier(Oswald) = 0.5
    curfac = max(0.9, 1.5*i1q/up)

    def label_line(ax, X, Y, U, V, label, color='k', size=8):
        """Add a label to a line, at the proper angle.

        Arguments
        ---------
        line : matplotlib.lines.Line2D object,
        label : str
        x : float
        x-position to place center of text (in data coordinated
        y : float
        y-position to place center of text (in data coordinates)
        color : str
        size : float
        """

        x1, x2 = X, X + U
        y1, y2 = Y, Y + V

        if y2 == 0:
            y2 = y1
        if x2 == 0:
            x2 = x1

        x = (x1 + x2) / 2
        y = (y1 + y2) / 2

        slope_degrees = np.rad2deg(np.angle(U + V * 1j))
        if slope_degrees < 0:
            slope_degrees += 180
        if 90 < slope_degrees <= 270:
            slope_degrees += 180

        x_offset = np.sin(np.deg2rad(slope_degrees))
        y_offset = np.cos(np.deg2rad(slope_degrees))
        bbox_props = dict(boxstyle="Round4, pad=0.1", fc="white", lw=0)
        text = ax.annotate(label, xy=(x, y), xytext=(x_offset * 10, y_offset * 8),
                           textcoords='offset points',
                           size=size, color=color,
                           horizontalalignment='center',
                           verticalalignment='center',
                           fontfamily='monospace', fontweight='bold', bbox=bbox_props)

        text.set_rotation(slope_degrees)
        return text

    if ax == 0:
        ax = plt.gca()
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    # ax.set_aspect('equal')

    ax.set_title(
        r'$U_1$={0} V, $I_1$={1} A, $U_p$={2} V'.format(
            round(u1, 1), round(i1, 1), round(up, 1)), fontsize=14)

    up /= uref
    ax.quiver(0, 0, 0, up, angles='xy', scale_units='xy', scale=1, units='dots',
              headwidth=qhw/2, headlength=qhl/2, headaxislength=qhl/2, width=qlw*2, color='k')
    label_line(ax, 0, 0, 0, up, '$U_p$', 'k', qts)

    ax.quiver(0, 0, u1d, u1q, angles='xy', scale_units='xy', scale=1, units='dots',
              headwidth=qhw, headlength=qhl, headaxislength=qhl, width=qlw, color='r')
    label_line(ax, 0, 0, u1d, u1q, '$U_1$', 'r', qts)

    ax.quiver(0, 1, uxd, 0, angles='xy', scale_units='xy', scale=1, units='dots',
              headwidth=qhw, headlength=qhl, headaxislength=qhl, width=qlw, color='g')
    label_line(ax, 0, 1, uxd, 0, '$U_d$', 'g', qts)

    ax.quiver(uxd, 1, 0, uxq, angles='xy', scale_units='xy', scale=1, units='dots',
              headwidth=qhw, headlength=qhl, headaxislength=qhl, width=qlw, color='g')
    label_line(ax, uxd, 1, 0, uxq, '$U_q$', 'g', qts)

    ax.quiver(0, 0, curfac*i1d, curfac*i1q, angles='xy', scale_units='xy', scale=1,
              units='dots', headwidth=qhw, headlength=qhl, headaxislength=qhl, width=qlw, color='b')
    label_line(ax, 0, 0, curfac*i1d, curfac*i1q, '$I_1$', 'b', qts)

    xmin, xmax = (min(0, uxd, i1d), max(0, i1d, uxd))
    ymin, ymax = (min(0, i1q, 1-uxq), max(1, i1q, 1+uxq))

    ax.set_xlim([xmin-0.1, xmax+0.1])
    ax.set_ylim([ymin-0.1, ymax+0.1])
    ax.grid(True)


def i1beta_phasor(up, i1, beta, r1, xd, xq, ax=0):
    """creates a phasor plot
    up: internal voltage
    i1: current
    beta: angle i1 vs up [deg]
    r1: resistance
    xd: reactance in direct axis
    xq: reactance in quadrature axis"""

    i1d, i1q = (i1*np.sin(beta/180*np.pi), i1*np.cos(beta/180*np.pi))
    uxdq = ((r1*i1d - xq*i1q), (r1*i1q + xd*i1d))
    __phasor_plot(ax, up, (i1d, i1q), uxdq)


def iqd_phasor(up, iqd, uqd, ax=0):
    """creates a phasor plot
    up: internal voltage
    iqd: current
    uqd: terminal voltage"""

    uxdq = (uqd[1]/np.sqrt(2), (uqd[0]/np.sqrt(2)-up))
    __phasor_plot(ax, up, (iqd[1]/np.sqrt(2), iqd[0]/np.sqrt(2)), uxdq)


def phasor(bch, ax=0):
    """create phasor plot from bch"""
    f1 = bch.machine['p']*bch.dqPar['speed']
    w1 = 2*np.pi*f1
    xd = w1*bch.dqPar['ld'][-1]
    xq = w1*bch.dqPar['lq'][-1]
    r1 = bch.machine['r1']
    if np.isscalar(bch.dqPar['up']):
        up = bch.dqPar['up']
    else:
        up = bch.dqPar['up'][-1]

    i1beta_phasor(up,
                  bch.dqPar['i1'][-1], bch.dqPar['beta'][-1],
                  r1, xd, xq, ax)
