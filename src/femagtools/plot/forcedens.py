"""
    femagtools.plot.forcedens
    ~~~~~~~~~~~~~~~~~~~~~~~~~

    Creating force density plots


"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

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


def forcedens(title, pos, fdens, ax=0):
    """plot force densities"""
    if ax == 0:
        ax = plt.gca()
    ax.set_title(title)
    ax.grid(True)

    ax.plot(pos, [1e-3*ft for ft in fdens[0]], label='F tang')
    ax.plot(pos, [1e-3*fn for fn in fdens[1]], label='F norm')
    ax.legend()
    ax.set_xlabel('Pos / deg')
    ax.set_ylabel('Force Density / kN/m²')


def forcedens_surface(fdens, ax=0):
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    xpos = [p for p in fdens.positions[0]['X']]
    ypos = [p['position'] for p in fdens.positions]
    z = 1e-3*np.array([p['FN']
                       for p in fdens.positions])
    _plot_surface(ax, xpos, ypos, z,
                  ('Rotor pos/°', 'Pos/°', 'F N / kN/m²'))


def forcedens_fft(title, fdens, nharm=40, ax=0):
    """plot force densities FFT
    Args:
      title: plot title
      fdens: force density object
      nharm: (int) num harmonics
    """
    if ax == 0:
        ax = plt.axes(projection="3d")

    F = 1e-3*fdens.fft(nharm)['fn_harm']['amplitude']
    fmin = 0.2
    num_bars = F.shape[0] + 1
    _xx, _yy = np.meshgrid(np.arange(1, num_bars),
                           np.arange(1, num_bars))
    z_size = F[F > fmin]
    x_pos, y_pos = _xx[F > fmin], _yy[F > fmin]
    z_pos = np.zeros_like(z_size)
    x_size = 2
    y_size = 2

    ax.bar3d(x_pos, y_pos, z_pos, x_size, y_size, z_size)
    #ax.view_init(azim=120)
    ax.set_xlim(0, num_bars+1)
    ax.set_ylim(0, num_bars+1)
    ax.set_title(title)
    ax.set_xlabel('M')
    ax.set_ylabel('N')
    ax.set_zlabel('kN/m²')
