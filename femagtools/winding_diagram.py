
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D


def _winding_data(Q, p, m):
    slots_per_pole_and_phase = int(Q / p / 2 / m)

    keys = [p + 1 for p in range(m) for _ in range(slots_per_pole_and_phase)] * p * 2

    direction = -1
    for i in range(Q):
        if i % slots_per_pole_and_phase == 0:
            direction *= -1
        keys[i] *= direction

    return keys


def winding_diagram(Q, p, m, filename):

    def color(key):
        colors = {
            1: "lime",
            2: "magenta",
            3: "gold",
        }
        try:
            return colors[key]
        except KeyError:
            return "blue"

    coil_height = 2.5
    thingy_height = 0.5
    thingy2_height = 0.3
    base_gap = 0.3
    arrow_head_length = .2

    data = _winding_data(Q, p, m)

    fig, ax = plt.subplots()

    for i, key in enumerate(data):

        coil_pos = i
        coilspan = int(Q / p / 2)
        base = -(abs(key) * base_gap + thingy_height + thingy2_height)

        out = key > 0
        coil_color = color(abs(key))

        if out:
            xdata = [
                coil_pos + coilspan / 2 - 0.1,
                coil_pos,
                coil_pos,
                coil_pos + coilspan / 2,
            ]
            ydata = [
                -thingy_height,
                0,
                coil_height,
                coil_height + thingy_height
            ]
            up_or_down = arrow_head_length
            arrow_y_pos = coil_height * .8
        else:
            xdata = [
                coil_pos - coilspan / 2,
                coil_pos,
                coil_pos,
                coil_pos - coilspan / 2 + 0.1,
            ]
            ydata = [
                coil_height + thingy_height,
                coil_height,
                0,
                -thingy_height,
            ]
            up_or_down = -arrow_head_length
            arrow_y_pos = coil_height * .2

        ax.arrow(coil_pos, arrow_y_pos, 0, up_or_down,
                 length_includes_head=True,
                 head_starts_at_zero=False,
                 head_length=arrow_head_length,
                 head_width=.2,
                 fc=coil_color,
                 lw=0,
                 )

        ax.add_line(Line2D(xdata, ydata, color=coil_color, lw=.8))
        ax.text(coil_pos, coil_height / 2,
                str(coil_pos + 1),
                fontsize=6,
                horizontalalignment="center",
                verticalalignment="center",
                backgroundcolor="white",
                bbox=dict(boxstyle='circle,pad=0', fc="white", lw=0),
                )

        ax.add_patch(Rectangle((coil_pos + 0.25, 0.1), 0.5, coil_height - 0.2, fc="lightblue"))

        if i == data.index(abs(key)) and abs(key) != 2:
            if out:
                x = coil_pos + coilspan / 2 - 0.1
                name = str(abs(key))
            else:
                x = coil_pos - coilspan / 2 + 0.1
                name = str(abs(key)) + "'"

            ax.add_line(Line2D([x, x], [-thingy_height, -2], color=coil_color, lw=.8))
            ax.text(x + 0.1, -2.2, name, color=coil_color)

        elif i == len(data) - 1 - data[::-1].index(-abs(key)) and abs(key) != 2:
            if out:
                x = coil_pos + coilspan / 2 - 0.1
            else:
                x = coil_pos - coilspan / 2 + 0.1

            ax.add_line(Line2D([x, x], [-thingy_height, -2], color=coil_color, lw=.8))
            ax.text(x + 0.1, -2.2, str(abs(key)) + "'", color=coil_color)

        elif abs(key) == 2 and (i == data.index(key) or data[data.index(key):].index(key)):
            if out:
                x = coil_pos + coilspan / 2 - 0.1
                name = str(abs(key))
            else:
                x = coil_pos - coilspan / 2 + 0.1
                name = str(abs(key)) + "'"

            ax.add_line(Line2D([x, x], [-thingy_height, -2], color=coil_color, lw=.8))
            ax.text(x + 0.1, -2.2, name, color=coil_color)

        else:
            if out:
                xdata = [
                    coil_pos + coilspan / 2 - 0.1,
                    coil_pos + coilspan / 2 - 0.1,
                    coil_pos - coilspan / 2,
                ]
            else:
                xdata = [
                    coil_pos - coilspan / 2 + 0.1,
                    coil_pos - coilspan / 2 + 0.1,
                    coil_pos + coilspan / 2,
                ]

            ydata = [
                -thingy_height,
                base,
                base
            ]
            ax.add_line(Line2D(xdata, ydata, color=coil_color, lw=.8))

    ax.autoscale(enable=True)
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.tight_layout()

    plt.savefig(filename)
