"""plot nc model

"""
# https://stackoverflow.com/questions/42386372/increase-the-speed-of-redrawing-contour-plot-in-matplotlib/42398244#42398244
import numpy as np
import matplotlib.pyplot as plt
import logging


logger = logging.getLogger("femagtools.plot.nc")


DEFAULT_CMAP='viridis'
"""default colormap (see https://matplotlib.org/stable/users/explain/colors/colormaps.html)"""


def spel(isa, superelements=[], with_axis=False, with_wiredir=False, ax=0):
    """plot super elements of I7/ISA7 model
    Args:
      isa: Isa7 object
      superelements: list of super elements (all if empty)
    """
    from matplotlib.patches import Polygon
    if ax == 0:
        ax = plt.gca()
    ax.set_aspect('equal')
    if superelements:
        spels = superelements
    else:
        spels = isa.superelements
    for se in spels:
        ax.add_patch(Polygon([n.xy
                            for nc in se.nodechains
                            for n in nc.nodes],
                            color=isa.color[se.color], lw=0))
        try:
            # draw wire direction
            if se.subregion and with_wiredir:
                if se.subregion.curdir != 0:
                        wkey = se.subregion.winding.key
                        if se.subregion.curdir < 0:
                            label = str(se.subregion.curdir*wkey)
                        else:
                            label = '+'+str(se.subregion.curdir*wkey)

                        xy = np.array([n.xy
                            for nc in se.nodechains
                            for n in nc.nodes])
                        cx, cy = np.mean(np.unique(xy[:, 0])), np.mean(np.unique(xy[:, -1]))
                        ax.text(cx, cy, label, rotation=np.arctan2(cy, cx)/np.pi*180-90,
                                horizontalalignment='center', verticalalignment='center')
        except:
            pass
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


def _contour(ax, title, elements, values, label='',
             cmap=DEFAULT_CMAP, isa=None, alpha=1):
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    if ax == 0:
        ax = plt.gca()
    ax.set_aspect('equal')
    if title:
        ax.set_title(title)
    if isa:
        for se in isa.superelements:
            ax.add_patch(Polygon([n.xy
                                  for nc in se.nodechains
                                  for n in nc.nodes],
                                 color='gray', alpha=0.1, lw=0))
    valid_values = np.logical_not(np.isnan(values))
    vertices = [[v.xy for v in e.vertices] for e in elements]
    patches = np.array([Polygon(xy) for xy in vertices])[valid_values]
    p = PatchCollection(patches, match_original=False,
                        cmap=cmap, alpha=alpha)
    p.set_array(np.asarray(values)[valid_values])
    ax.add_collection(p)
    cb = plt.colorbar(p, shrink=0.9)

    for patch in np.array([Polygon(xy, fc='white', alpha=1.0)
                           for xy in vertices])[np.isnan(values)]:
        ax.add_patch(patch)
    if label:
        cb.set_label(label=label)
    ax.autoscale(enable=True)
    ax.axis('off')


def demag(isa, cmap=DEFAULT_CMAP, ax=0):
    """plot demag of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    emag = [e for e in isa.elements if e.is_magnet()]
    demag = np.array([e.demagnetization(isa.MAGN_TEMPERATURE) for e in emag])
    _contour(ax, f'Demagnetization at {isa.MAGN_TEMPERATURE} °C (max -{np.max(demag):.1f} kA/m)',
             emag, demag, '-H / kA/m', cmap, isa)
    logger.info("Max demagnetization %f", np.max(demag))


def remanence(isa, cmap=DEFAULT_CMAP, ax=0):
    """plot remanence of NC/I7/ISA7 model
    Args:
      isa: Isa7/NC object
    """
    emag = [e for e in isa.elements if e.is_magnet()]
    rem = np.linalg.norm([e.remanence(isa.MAGN_TEMPERATURE)
                          for e in emag], axis=1)
    _contour(ax, f'Remanence at {isa.MAGN_TEMPERATURE} °C (min {np.min(rem):.1f} T)',
             emag, rem, 'T', cmap, isa)
    logger.info("Min remanence %f", np.min(rem))


def demag_pos(isa, pos=-1, icur=-1, ibeta=-1, cmap=DEFAULT_CMAP, ax=0):
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
    hmax = np.max(hpol)
    hpol[hpol == 0] = np.nan
    _contour(ax, f'Demagnetization at pos. {round(x/np.pi*180):.1f}°, '
    f'{isa.MAGN_TEMPERATURE} °C (max -{hmax:.1f} kA/m)',
             emag, hpol, '-H / kA/m', cmap, isa)
    logger.info("Max demagnetization %f kA/m", np.nanmax(hpol))


def __elements_of_subreg(isa, subreg):
    if subreg:
        if isinstance(subreg, list):
            sr = subreg
        else:
            sr = [subreg]
        for s in sr:
            yield from isa.get_subregion(s).elements()
    else:
        for e in isa.elements:
            yield e


def flux_density(isa, subreg=[], cmap=DEFAULT_CMAP, ax=0):
    """plot flux density of NC/I7/ISA7 model

    Args:
        isa: Isa7/NC object
        subreg: list of subregion names (all if empty)
        icur: cur index or last index if -1
        ibeta: beta angle index or last index if -1
    """
    elements = [e for e in __elements_of_subreg(isa, subreg)]
    fluxd = np.array([np.linalg.norm(e.flux_density()) for e in elements])
    _contour(ax, f'Flux Density T (max {np.max(fluxd):.1f} T)',
             elements, fluxd, '', cmap)
    logger.info("Max flux dens %f", np.max(fluxd))


def flux_density_eccentricity(isa, subreg=[], icur=-1, ibeta=-1,
                              cmap='plasma_r', ax=0, alpha=0.75):
    """plot eccentricity (axis ratio) of flux density in lamination

    Args:
        isa: Isa7/NC object
        subreg: list of subregion names (all if empty)
        icur: cur amplitude index or last index if -1
        ibeta: beta angle index or last index if -1
    """
    from ..utils import fft
    elements = []
    ecc = []
    pulsating = 0
    # eliminate double values at end
    # TODO: adapt to linear machines
    pos = [p
           for p in isa.pos_el_fe_induction
           if p < 2*np.pi/isa.pole_pairs] + [2*np.pi/isa.pole_pairs]
    i = len(pos)
    apos = np.array(pos)/np.pi*180
    for e in __elements_of_subreg(isa, subreg):
        if e.is_lamination():
            br = isa.el_fe_induction_1[e.key-1, :i+1, icur, ibeta]
            bt = isa.el_fe_induction_2[e.key-1, :i+1, icur, ibeta]
            brtmax = np.max(br-np.mean(br)), np.max(bt-np.mean(bt))
            if np.all(np.isclose(brtmax, 0)):
                continue
            elements.append(e)
            if np.any(np.isclose(brtmax, 0)):
                ecc.append(pulsating)
            else:
                br0 = fft(apos, br-np.mean(br))
                x = br0['a']*np.cos(2*np.pi*apos/br0['T0']+br0['alfa0'])
                bt0 = fft(apos, bt-np.mean(bt))
                y = bt0['a']*np.cos(2*np.pi*apos/bt0['T0']+bt0['alfa0'])
                if (br0['a'] > br0['nue'][isa.pole_pairs]
                    or bt0['a'] > bt0['nue'][isa.pole_pairs]):
                    ecc.append(pulsating)
                else:
                    kmax = np.argmax(np.linalg.norm((x, y), axis=0))
                    kmin = np.argmin(np.linalg.norm((x, y), axis=0))
                    a = np.linalg.norm((x[kmax], y[kmax]))
                    b = np.linalg.norm((x[kmin], y[kmin]))
                    ecc.append(b/a) #np.sqrt(1-b**2/a**2))

    _contour(ax, '', #'Eccentricity of Flux Density',
                 elements, ecc, 'axis ratio', cmap, alpha=alpha)


def flux_density_pos(isa, ipos, subreg=[], icur=-1, ibeta=-1, cmap=DEFAULT_CMAP, ax=0):
    """plot flux density at rotor pos for each element of NC/I7/ISA7 model

    Args:
        isa: Isa7/NC object
        ipos: position index
        icur: cur index or last index if -1
        ibeta: beta angle index or last index if -1

    """
    elements = [e for e in __elements_of_subreg(isa, subreg)
                if e not in isa.airgap_center_elements]
    b = []
    for e in elements:
        fd = isa.flux_density(e, icur, ibeta)
        b.append(np.linalg.norm(
            (fd['bx'][ipos], fd['by'][ipos])))
    fluxd = np.array(b)
    pos = isa.pos_el_fe_induction[ipos]*180/np.pi
    isa.rotate(isa.pos_el_fe_induction[ipos])
    _contour(ax, f'Flux Density T at {pos:.1f}° (max {np.max(fluxd):.1f} T)',
             elements, fluxd, '', cmap)
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


def loss_density(isa, subreg=[], cmap=DEFAULT_CMAP, ax=0):
    """plot loss density of NC/I7/ISA7 model

    Args:
        isa: Isa7/NC object
    """
    elements = [e for e in __elements_of_subreg(isa, subreg)]
    lossd = np.array([e.loss_density*1e-3 for e in elements])
    _contour(ax, 'Loss Density kW/m³', elements, lossd, '', cmap)

def temperature_distribution(isa, ax=0, cmap='plasma'):
    """plot temperature distribution of NC/I7/ISA7 model

    Args:
        isa: Isa7/NC object
    """
    temp = []
    elements = [e for e in isa.elements]
    for e in isa.elements:
        tmp = 0
        ctr = 1
        for n in e.vertices:
            tmp += n.vpot[-1]
            ctr = ctr + 1
        temp.append(tmp/ctr)
    _contour(ax, 'Temperature in K', elements, temp, '', cmap)
