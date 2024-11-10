"""
    femagtools.plot.mcv
    ~~~~~~~~~~~~~~~~~~~

    Creating mcv plots


"""
import numpy as np
import matplotlib.pyplot as plt


def mcv_hbj(mcv, log=True, ax=0):
    """plot H, B, J of mcv dict"""
    import femagtools.mcv
    MUE0 = 4e-7*np.pi
    ji = []

    csiz = len(mcv['curve'])
    if ax == 0:
        ax = plt.gca()
    ax.set_title(mcv['name'])
    for k, c in enumerate(mcv['curve']):
        bh = [(bi, hi*1e-3)
              for bi, hi in zip(c['bi'],
                                c['hi'])]
        try:
            if csiz == 1 and mcv['ctype'] in (femagtools.mcv.MAGCRV,
                                              femagtools.mcv.ORIENT_CRV):
                ji = [b-MUE0*h*1e3 for b, h in bh]
        except Exception:
            pass
        bi, hi = zip(*bh)

        label = 'Flux Density'
        if csiz > 1:
            label = 'Flux Density ({0}°)'.format(mcv.mc1_angle[k])
        if log:
            ax.semilogx(hi, bi, label=label)
            ax.semilogx(hi, bi, 'k.')
            if ji:
                ax.semilogx(hi, ji, label='Polarisation')
        else:
            ax.plot(hi, bi, label=label)
            ax.plot(hi, bi, 'k.')
            if ji:
                ax.plot(hi, ji, label='Polarisation')
    ax.set_xlabel('H / kA/m')
    ax.set_ylabel('T')
    if ji or csiz > 1:
        ax.legend(loc='lower right')
    ax.grid()


def mcv_muer(mcv, ax=0):
    """plot rel. permeability vs. B of mcv dict"""
    MUE0 = 4e-7*np.pi
    bi, ur = zip(*[(bx, bx/hx/MUE0)
                   for bx, hx in zip(mcv['curve'][0]['bi'],
                                     mcv['curve'][0]['hi']) if not hx == 0])
    if ax == 0:
        ax = plt.gca()
    ax.plot(bi, ur)
    ax.plot(bi, ur, 'k.')
    ax.set_xlabel('B / T')
    ax.set_title('rel. Permeability')
    ax.grid()

def mcv_nuer(mcv, ax=0):
    """plot rel. reluctivity vs. B of mcv dict"""
    MUE0 = 4e-7*np.pi
    bi, ur = zip(*[(bx, bx/hx/MUE0)
                   for bx, hx in zip(mcv['curve'][0]['bi'],
                                     mcv['curve'][0]['hi']) if not hx == 0])
    if ax == 0:
        ax = plt.gca()
    nuer = 1/np.array(ur)
    ax.plot(bi, 1/np.array(ur))
    ax.plot(bi, nuer, 'k.')
    ax.set_xlabel('B / T')
    ax.set_title('rel. Reluctivity')
    ax.grid()


def felosses(losses, coeffs, title='', log=True, ax=0):
    """plot iron losses with steinmetz or jordan approximation

    Args:
      losses: dict with f, B, pfe values
      coeffs: list with steinmetz (cw, alpha, beta) or
              jordan (cw, alpha, ch, beta, gamma)  or
              bertotti (cw, alpha, cw, ce) coeffs
      title: title string
      log: log scale for x and y axes if True

    """
    import femagtools.losscoeffs as lc
    if ax == 0:
        ax = plt.gca()

    fo = losses['fo']
    Bo = losses['Bo']
    B = np.linspace(0.9*np.min(losses['B']),
                        1.1*0.9*np.max(losses['B']))

    for i, f in enumerate(losses['f']):
        pfe = [p for p in np.array(losses['pfe'])[i] if p]
        if f > 0:
            if len(coeffs) == 5:
                ax.plot(B, lc.pfe_jordan(f, B, *coeffs, fo=fo, Bo=Bo))
            if len(coeffs) == 4:
                ax.plot(B, lc.pfe_bertotti(f, B, *coeffs))
            elif len(coeffs) == 3:
                ax.plot(B, lc.pfe_steinmetz(f, B, *coeffs, fo=fo, Bo=Bo))
            ax.plot(losses['B'][:len(pfe)], pfe,
                    marker='o', label="{} Hz".format(f))

    ax.set_title("Fe Losses/(W/kg) " + title)
    if log:
        ax.set_yscale('log')
        ax.set_xscale('log')
    ax.set_xlabel("Flux Density [T]")
    # plt.ylabel("Pfe [W/kg]")
    ax.legend()
    ax.grid(True)
