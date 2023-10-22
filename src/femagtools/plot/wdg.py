"""
    femagtools.plot.wdg
    ~~~~~~~~~~~~~~~~~~~

    Creating winding plots

"""
import numpy as np
import matplotlib.pyplot as plt


def mmf(f, title='', ax=0):
    """plot magnetomotive force (mmf) of winding"""
    if ax == 0:
        ax = plt.gca()
    if title:
        ax.set_title(title)
    ax.plot(np.array(f['pos'])/np.pi*180, f['mmf'])
    ax.plot(np.array(f['pos_fft'])/np.pi*180, f['mmf_fft'])
    ax.set_xlabel('Position / Deg')

    phi = [f['alfa0']/np.pi*180, f['alfa0']/np.pi*180]
    y = [min(f['mmf_fft']), 1.1*max(f['mmf_fft'])]
    ax.plot(phi, y, '--')
    alfa0 = round(f['alfa0']/np.pi*180, 2)
    ax.text(phi[0]/2, y[0]+0.05, f"{alfa0}°",
            ha="center", va="bottom")
    ax.annotate(f"", xy=(phi[0], y[0]),
                xytext=(0, y[0]), arrowprops=dict(arrowstyle="->"))
    ax.grid()


def mmf_fft(f, title='', mmfmin=1e-2, ax=0):
    """plot winding mmf harmonics"""
    if ax == 0:
        ax = plt.gca()
    if title:
        ax.set_title(title)
    else:
        ax.set_title('MMF Harmonics (per phase)')
    ax.grid(True)
    order, mmf = np.array([(n, m) for n, m in zip(f['nue'],
                                                  f['mmf_nue']) if m > mmfmin]).T
    try:
        markerline1, stemlines1, _ = ax.stem(order, mmf, '-.', basefmt=" ")
        ax.set_xticks(order)
    except ValueError:  # empty sequence
        pass


def zoneplan(wdg, ax=0):
    """plot zone plan of winding wdg"""
    from matplotlib.patches import Rectangle
    upper, lower = wdg.zoneplan()
    Qb = len([n for l in upper for n in l])
    from femagtools.windings import coil_color
    rh = 0.5
    if lower:
        yl = rh
        ymax = 2*rh + 0.2
    else:
        yl = 0
        ymax = rh + 0.2
    if ax == 0:
        ax = plt.gca()
    ax.axis('off')
    ax.set_xlim([-0.5, Qb-0.5])
    ax.set_ylim([0, ymax])
    ax.set_aspect(Qb/6+0.3)

    for i, p in enumerate(upper):
        for x in p:
            ax.add_patch(Rectangle((abs(x)-1.5, yl), 1, rh,
                                   facecolor=coil_color[i],
                                   edgecolor='white', fill=True))
            s = f'+{i+1}' if x > 0 else f'-{i+1}'
            ax.text(abs(x)-1, yl+rh/2, s, color='black',
                    ha="center", va="center")
    for i, p in enumerate(lower):
        for x in p:
            ax.add_patch(Rectangle((abs(x)-1.5, yl-rh), 1, rh,
                                   facecolor=coil_color[i],
                                   edgecolor='white', fill=True))
            s = f'+{i+1}' if x > 0 else f'-{i+1}'
            ax.text(abs(x)-1, yl-rh/2, s, color='black',
                    ha="center", va="center")

    yu = yl+rh
    step = 1 if Qb < 25 else 2
    if lower:
        yl -= rh
    margin = 0.05
    ax.text(-0.5, yu+margin, f'Q={wdg.Q}, p={wdg.p}, q={round(wdg.q,4)}',
            ha='left', va='bottom', size=15)
    for i in range(0, Qb, step):
        ax.text(i, yl-margin, f'{i+1}', ha="center", va="top")


def winding_factors(wdg, n=8, ax=0):
    """plot winding factors"""
    ax = plt.gca()
    ax.set_title(f'Winding factors Q={wdg.Q}, p={wdg.p}, q={round(wdg.q,4)}')
    ax.grid(True)
    order, kwp, kwd, kw = np.array(
        [(n, k1, k2, k3)
         for n, k1, k2, k3 in zip(wdg.kw_order(n),
                                  wdg.kwp(n),
                                  wdg.kwd(n),
                                  wdg.kw(n))]).T
    try:
        markerline1, stemlines1, _ = ax.stem(order-1, kwp,
                                             'C1:', basefmt=" ",
                                             markerfmt='C1.',
                                             label='Pitch')
        markerline2, stemlines2, _ = ax.stem(order+1, kwd,
                                             'C2:', basefmt=" ",
                                             markerfmt='C2.',
                                             label='Distribution')
        markerline3, stemlines3, _ = ax.stem(order, kw,
                                             'C0-', basefmt=" ",
                                             markerfmt='C0o',
                                             label='Total')
        ax.set_xticks(order)
        ax.legend()
    except ValueError:  # empty sequence
        pass


def winding(wdg, ax=0):
    """plot coils of windings wdg"""
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from femagtools.windings import coil_color

    coil_len = 25
    coil_height = 4
    dslot = 8
    arrow_head_length = 2
    arrow_head_width = 2

    if ax == 0:
        ax = plt.gca()
    z = wdg.zoneplan()
    xoff = 0
    if z[-1]:
        xoff = 0.75
    yd = dslot*wdg.yd
    mh = 2*coil_height/yd
    slots = sorted([abs(n) for m in z[0] for n in m])
    smax = slots[-1]*dslot
    for n in slots:  # draw slots and lamination
        x = n*dslot
        ax.add_patch(Rectangle((x + dslot/4, 1), dslot /
                     2, coil_len - 2, fc="lightblue"))
        ax.text(x, coil_len / 2,
                str(n),
                horizontalalignment="center",
                verticalalignment="center",
                backgroundcolor="white",
                bbox=dict(boxstyle='circle,pad=0', fc="white", lw=0))

    nl = 2 if z[1] else 1
    line_thickness = [0.6, 1.2]
    for i, layer in enumerate(z):
        b = -xoff if i else xoff
        lw = line_thickness[i]
        direction = ['right', 'left']
        d = 1
        for m, mslots in enumerate(layer):
            for k in mslots:
                x = abs(k) * dslot + b
                kcoil = (abs(k)-1)//wdg.yd
                xpoints = []
                ypoints = []
                if nl == 2:
                    if k > 0:
                        d = 0 if i == 0 else 1
                    else:
                        d = 1 if i == 1 else 0
                else:
                    d = 0 if k > 0 else 1

                if direction[d] == 'right':
                    # first layer, positive dir or neg. dir and 2-layers:
                    #   from right bottom
                    if x + yd > smax+b:
                        dx = dslot if yd > dslot else yd/4
                        xpoints = [x + yd//2 + dx - xoff]
                        ypoints = [-coil_height + mh*dx]
                    xpoints += [x + yd//2 - xoff, x, x, x + yd//2-xoff]
                    ypoints += [-coil_height, 0, coil_len,
                                coil_len+coil_height]
                    if x + yd > smax+b:
                        xpoints += [x + yd//2 + dx - xoff]
                        ypoints += [coil_len+coil_height - mh*dx]
                else:
                    # from left bottom
                    if x - yd < 0:  # and x - yd/2 > -3*dslot:
                        dx = dslot if yd > dslot else yd/4
                        xpoints = [x - yd//2 - dx + xoff]
                        ypoints = [- coil_height + mh*dx]
                    xpoints += [x - yd//2+xoff, x, x, x - yd/2+xoff]
                    ypoints += [-coil_height, 0, coil_len,
                                coil_len+coil_height]
                    if x - yd < 0:  # and x - yd > -3*dslot:
                        xpoints += [x - yd//2 - dx + xoff]
                        ypoints += [coil_len + coil_height - mh*dx]

                ax.add_line(Line2D(xpoints, ypoints,
                            color=coil_color[m], lw=lw))

                if k > 0:
                    h = arrow_head_length
                    y = coil_len * 0.8
                else:
                    h = -arrow_head_length
                    y = coil_len * 0.2
                ax.arrow(x, y, 0, h,
                         length_includes_head=True,
                         head_starts_at_zero=False,
                         head_length=arrow_head_length,
                         head_width=arrow_head_width,
                         fc=coil_color[m], lw=0)
    if False:  # TODO show winding connections
        m = 0
        for k in [n*wdg.Q/wdg.p/wdg.m + 1 for n in range(wdg.m)]:
            if k < len(slots):
                x = k * dslot + b + yd/2 - xoff
                ax.add_line(Line2D([x, x],
                                   [-2*coil_height, -coil_height],
                                   color=coil_color[m], lw=lw))
                ax.text(x, -2*coil_height+0.5, str(m+1), color=coil_color[m])
            m += 1
    ax.autoscale(enable=True)
    ax.set_axis_off()
