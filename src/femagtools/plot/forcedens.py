"""
    femagtools.plot.forcedens
    ~~~~~~~~~~~~~~~~~~~~~~~~~

    Creating force density plots


"""
import numpy as np
import matplotlib.pyplot as plt
try:
    from matplotlib import colormaps
    cmap_viridis = colormaps['viridis']
    cmap_jet = colormaps['jet']
    cmap_YlOrBr = colormaps['YlOrBr']
except:  # older matplotlib
    from matplotlib import cm
    cmap_viridis = cm.viridis
    cmap_jet = cm.jet
    cmap_YlOrBr = cm.YlOrBr


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


def _plot_surface(ax, x, y, z, labels, azim=None, cmap=cmap_viridis):
    """helper function for surface plots"""
    # ax.tick_params(axis='both', which='major', pad=-3)
    assert np.size(x) > 1 and np.size(y) > 1 and np.size(z) > 1
    if azim is not None:
        ax.azim = azim
    X, Y = np.meshgrid(x, y)
    Z = np.ma.masked_invalid(z)
    ax.plot_surface(X, Y, Z,
                    rstride=1, cstride=1,
                    cmap=cmap,
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


def forcedens_surface(fdens, attr='FN', ax=0, cmap=cmap_jet):
    if ax == 0:
        _create_3d_axis()
        ax = plt.gca()
    z = 1e-3*np.array(getattr(fdens, attr))
    n, m = z.shape
    x = np.linspace(0, 360, m)
    y = np.linspace(0, 360, n)
    _plot_surface(ax, x, y, z,
                  ('Rotor pos/°', 'Pos/°', f'{attr} / kN/m²'),
                  cmap=cmap)

def forcedens_contour(fdens, attr='FN', ax=0, cmap=cmap_jet):
    if ax == 0:
        ax = plt.gca()
    z = 1e-3*np.array(getattr(fdens, attr))
    n, m = z.shape
    x = np.linspace(0, 360, m)
    y = np.linspace(0, 360, n)
    cs = ax.contourf(x, y, z, cmap=cmap)
    ax.set_xlabel('Rotor pos/°')
    ax.set_ylabel('Pos/°')
    if attr == 'FN':
        ax.set_title('Radial Force Density')
    elif attr == 'FT':
        ax.set_title('Tangential Force Density')
    cfig = ax.get_figure()
    cbar = cfig.colorbar(cs, ax=ax)
    cbar.ax.set_ylabel('kN/m²')

def forcedens_fft(title, fdens, harmmax=(), #(200, 40),
                  cmap=cmap_YlOrBr,
                  ax=0):
    """plot force densities FFT
    Args:
      title: plot title
      fdens: force density object
      cmap: colormap
    """
    if ax == 0:
        ax = plt.axes(projection="3d")

    FN = 1e-3*fdens.fft()['fn_harm']['amplitude']
    if harmmax:
        shape = np.min((harmmax, np.shape(FN)), axis=0)
    else:
        shape = np.shape(FN)

    _x = np.arange(shape[0])+1
    _y = np.arange(shape[1])+1
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    top = np.ravel(FN[:shape[0],:shape[1]])
    bottom = np.zeros_like(top)
    maxy = max(x[-1], y[-1])
    width, depth = 5*x[-1]/maxy, 5*y[-1]/maxy

    min_height = np.min(top)
    max_height = np.max(top)
    rgba = [cmap((k-min_height)/max_height) for k in top]
    ax.bar3d(x, y, bottom, width, depth, top, color=rgba)

    ax.set_xlim(0, shape[0]+1)
    ax.set_ylim(0, shape[1]+1)
    ax.set_title(title)
    ax.set_xlabel('M')
    ax.set_ylabel('N')
    ax.set_zlabel('kN/m²')
