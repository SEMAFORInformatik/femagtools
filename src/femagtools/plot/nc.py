"""plot nc model

"""
# https://stackoverflow.com/questions/42386372/increase-the-speed-of-redrawing-contour-plot-in-matplotlib/42398244#42398244
import numpy as np
import matplotlib.pyplot as plt
import logging


logger = logging.getLogger("femagtools.plot.nc")


def spel(isa, with_axis=False, ax=0):
    """plot super elements of I7/ISA7 model
    Args:
      isa: Isa7 object
    """
    from matplotlib.patches import Polygon
    if ax == 0:
        ax = plt.gca()
    ax.set_aspect('equal')
    for se in isa.superelements:
        ax.add_patch(Polygon([n.xy
                              for nc in se.nodechains
                              for n in nc.nodes],
                             color=isa.color[se.color], lw=0))

    ax.autoscale(enable=True)
    if not with_axis:
        ax.axis('off')


def mesh(isa, with_axis=False, ax=0):
    """plot mesh of I7/ISA7 model
    Args:
      isa: Isa7 object
    """
    from matplotlib.lines import Line2D
    if ax == 0:
        ax = plt.gca()
    ax.set_aspect('equal')
    for el in isa.elements:
        z = np.array([v.xy for v in el.vertices])
        pts = np.vstack((z, z[0])).T
        ax.add_line(Line2D(pts[0], pts[1],
                           color='b', ls='-', lw=0.25))

    # for nc in isa.nodechains:
    #    pts = [list(i) for i in zip(*[(n.x, n.y) for n in nc.nodes])]
    #    ax.add_line(Line2D(pts[0], pts[1], color="b", ls="-", lw=0.25,
    #                       marker=".", ms="2", mec="None"))

    # for nc in isa.nodechains:
    #    if nc.nodemid is not None:
    #        plt.plot(*nc.nodemid.xy, "rx")

    ax.autoscale(enable=True)
    if not with_axis:
        ax.axis('off')


def _contour(ax, title, elements, values, label='', isa=None):
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    if ax == 0:
        ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_title(title)
    if isa:
        for se in isa.superelements:
            ax.add_patch(Polygon([n.xy
                                  for nc in se.nodechains
                                  for n in nc.nodes],
                                 color='gray', alpha=0.1, lw=0))
    valid_values = np.logical_not(np.isnan(values))
    patches = np.array([Polygon([v.xy for v in e.vertices])
                       for e in elements])[valid_values]
    # , cmap=matplotlib.cm.jet, alpha=0.4)
    p = PatchCollection(patches, alpha=1.0, match_original=False)
    p.set_array(np.asarray(values)[valid_values])
    ax.add_collection(p)
    cb = plt.colorbar(p, shrink=0.9)
    for patch in np.array([Polygon([v.xy for v in e.vertices],
                                   fc='white', alpha=1.0)
                           for e in elements])[np.isnan(values)]:
        ax.add_patch(patch)
    if label:
        cb.set_label(label=label)
    ax.autoscale(enable=True)
    ax.axis('off')


def demag(isa, ax=0):
    """plot demag of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    emag = [e for e in isa.elements if e.is_magnet()]
    demag = np.array([e.demagnetization(isa.MAGN_TEMPERATURE) for e in emag])
    _contour(ax, f'Demagnetization at {isa.MAGN_TEMPERATURE} °C (max -{np.max(demag):.1f} kA/m)',
             emag, demag, '-H / kA/m', isa)
    logger.info("Max demagnetization %f", np.max(demag))


def demag_pos(isa, pos=-1, icur=-1, ibeta=-1, ax=0):
    """plot demag of NC/I7/ISA7 model at rotor position
    Args:
      isa: Isa7/NC object
      pos: rotor position in degree (maximum h if -1)
      icur: cur amplitude index or last index if -1
      ibeta: beta angle index or last index if -1
    """
    emag = [e for e in isa.elements if e.is_magnet()]
    demag = np.array([isa.demagnetization(e, icur, ibeta)[1]
                      for e in emag])
    if pos>=0:
        for i, x in enumerate(isa.pos_el_fe_induction):
            if x >= pos/180*np.pi:
                break
    else:
        demagmax = np.max(demag, axis=0)
        i = np.argmax(demagmax)
        x = isa.pos_el_fe_induction[i]

    hpol = demag[:, i]
    hpol[hpol == 0] = np.nan
    _contour(ax, f'Demagnetization at pos. {round(x/np.pi*180):.1f}°,'
    f'{isa.MAGN_TEMPERATURE} °C (max -{np.max(hpol):.1f} kA/m)',
             emag, hpol, '-H / kA/m', isa)
    logger.info("Max demagnetization %f kA/m", np.nanmax(hpol))


def __elements_of_subreg(isa, subreg):
    if subreg:
        if isinstance(subreg, list):
            sr = subreg
        else:
            sr = [subreg]
        for s in sr:
            for e in isa.get_subregion(s).elements():
                yield e
    else:
        for e in isa.elements:
            yield e


def flux_density(isa, subreg=[], ax=0):
    """plot flux density of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    elements = [e for e in __elements_of_subreg(isa, subreg)]
    fluxd = np.array([np.linalg.norm(e.flux_density()) for e in elements])
    _contour(ax, f'Flux Density T (max {np.max(fluxd):.1f} T)',
             elements, fluxd)
    logger.info("Max flux dens %f", np.max(fluxd))


def max_flux_density(isa, subreg=[], icur=-1, ibeta=-1, ax=0):
    """plot max flux density of each element of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    elements = [e for e in __elements_of_subreg(isa, subreg)]
    bmax = []
    for e in elements:
        fd = isa.flux_density(e, icur, ibeta)
        bmax.append(np.max(
            np.linalg.norm((fd['bx'], fd['by']), axis=0)))
    fluxd = np.array(bmax)
    _contour(ax, f'Max Flux Density T (max {np.max(fluxd):.1f} T)',
             elements, fluxd)
    logger.info("Max flux dens %f", np.max(fluxd))


def min_flux_density(isa, subreg=[], icur=-1, ibeta=-1, ax=0):
    """plot min flux density of each element of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    elements = [e for e in __elements_of_subreg(isa, subreg)]
    bmin = []
    for e in elements:
        fd = isa.flux_density(e, icur, ibeta)
        bmin.append(np.min(
            np.linalg.norm((fd['bx'], fd['by']), axis=0)))

    fluxd = np.array(bmin)
    _contour(ax, f'Min Flux Density T (max {np.max(fluxd):.1f} T)',
             elements, fluxd)
    logger.info("Max flux dens %f (element %d)",
                np.max(fluxd), elements[np.argmax(fluxd)].key)


def flux_density_pos(isa, ipos, subreg=[], icur=-1, ibeta=-1, ax=0):
    """plot flux density at rotor pos for each element of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    elements = [e for e in __elements_of_subreg(isa, subreg)]
    b = []
    for e in elements:
        fd = isa.flux_density(e, icur, ibeta)
        b.append(fd['bx'][ipos])
    fluxd = np.array(b)
    pos = isa.pos_el_fe_induction[ipos]*180/np.pi
    isa.rotate(isa.pos_el_fe_induction[ipos])
    _contour(ax, f'Flux Density T at {pos:.1f}° (max {np.max(fluxd):.1f} T)',
             elements, fluxd)
    logger.info("Max flux dens %f", np.max(fluxd))
    isa.rotate(0)


def airgap_flux_density_pos(isa, ipos, icur=-1, ibeta=-1, ax=0):
    """plot flux density at rotor pos for each element of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    bx = []
    by = []
    pos = []
    for e in isa.airgap_center_elements:
        fd = isa.flux_density(e, icur, ibeta)
        bx.append(fd['bx'][ipos])
        by.append(fd['by'][ipos])
        pos.append(np.arctan2(e.center[1], e.center[0])/np.pi*100)
    if ax == 0:
        ax = plt.gca()
    ax.plot(pos, bx)
    ax.plot(pos, by)
    ax.grid()

    logger.info("Max flux dens %f", np.max(np.abs(bx)))


def loss_density(isa, subreg=[], ax=0):
    """plot loss density of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    elements = [e for e in __elements_of_subreg(isa, subreg)]
    lossd = np.array([e.loss_density*1e-3 for e in elements])
    _contour(ax, 'Loss Density kW/m³', elements, lossd)
