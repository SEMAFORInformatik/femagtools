"""
Create longitudinal drawing of radial flux machine
"""
import matplotlib.pyplot as plt
import matplotlib.patches as pch

def _draw_shaft(ax, dy2, lfe):
    xx = (0.0, lfe, lfe, 0.0)
    yy = (dy2/2, dy2/2, -dy2/2, -dy2/2)
    ax.fill(xx, yy,
            facecolor='lightgrey',
            edgecolor='black',
            linewidth=0)
    xx = (-lfe/4, lfe+lfe/4, lfe+lfe/4, -lfe/4)
    yy = (dy2/4, dy2/4, -dy2/4, -dy2/4)
    ax.fill(xx, yy,
            facecolor='lightgrey',
            edgecolor='black',
            linewidth=0)

def _draw_rotor(ax, da1, dy2, lfe):
    ag = 0.02*da1
    xx = (0.0, lfe, lfe, 0.0)
    yy = (dy2/2, dy2/2, da1/2-ag, da1/2-ag)
    ax.fill(xx, yy,
            facecolor='skyblue',
            edgecolor='black',
            linewidth=0)
    yy = (-dy2/2, -dy2/2, -da1/2+ag, -da1/2+ag)
    ax.fill(xx, yy,
            facecolor='skyblue',
            edgecolor='black',
            linewidth=0)

def _draw_stator(ax, da1, dy1, lfe):
    xx = (0.0, lfe, lfe, 0.0)
    yy = (da1/2, da1/2, dy1/2, dy1/2)
    # yoke
    ax.fill(xx, yy,
            facecolor='skyblue',
            edgecolor='black',
            linewidth=0)

    yy = (-da1/2, -da1/2, -dy1/2, -dy1/2)
    ax.fill(xx, yy,
            facecolor='skyblue',
            edgecolor='black',
            linewidth=0)

    # winding
    yh = (dy1-da1)/2
    xx = (-yh/2, 0, 0, -yh/2)
    yy = (da1/2, da1/2, dy1/2-yh/2, dy1/2-yh/2)
    ax.fill(xx, yy, facecolor='gold',
                edgecolor='black',
                linewidth=0)

    xx = (lfe, lfe+yh/2, lfe+yh/2, lfe)
    ax.fill(xx, yy, facecolor='gold',
                edgecolor='black',
                linewidth=0)

    yy = (-da1/2, -da1/2, -dy1/2+yh/2, -dy1/2+yh/2)
    ax.fill(xx, yy, facecolor='gold',
                edgecolor='black',
                linewidth=0)

    xx = (-yh/2, 0, 0, -yh/2)
    ax.fill(xx, yy, facecolor='gold',
                edgecolor='black',
                linewidth=0)

def machine(machine, ax):
    dy2 = machine['inner_diam']*1e3
    dy1 = machine['outer_diam']*1e3
    da1 = machine['bore_diam']*1e3
    lfe = machine['lfe']*1e3

    _draw_rotor(ax, da1, dy2, lfe)
    _draw_stator(ax, da1, dy1, lfe)
    _draw_shaft(ax, dy2, lfe)
    ax.set_aspect('equal')

    for loc, spine in ax.spines.items():
        spine.set_color('none')  # don't draw spine
    #ax.yaxis.set_ticks([])
    #ax.xaxis.set_ticks([])


if __name__ == '__main__':
    machine1 = {
        "outer_diam": 0.2442,
        "bore_diam": 0.179,
        "inner_diam": 0.06,
        "airgap": 0.7e-3,
        "lfe": 0.083,
    }
    fig, ax = plt.subplots()
    machine(machine1, ax)
    plt.show()
