
from xml.etree import ElementTree as ET


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

    coilspan = int(Q / p / 2) * 10
    coil_height = 25
    top_height = 5
    neck_height = 3
    base_gap = 3
    arrow_head_length = 2
    arrow_head_width = 2

    width = Q * 10 + coilspan / 4 + 10
    min_x = -(coilspan / 4)

    svg = ET.Element("svg", dict(xmlns="http://www.w3.org/2000/svg", viewBox=f"{min_x} -40 {width} 70"))

    data = _winding_data(Q, p, m)

    for i, key in enumerate(data):

        coil_pos = i * 10
        base = abs(key) * base_gap + top_height + neck_height

        out = key > 0
        coil_color = color(abs(key))

        if out:
            xdata = [
                coil_pos + coilspan / 2 - 1,
                coil_pos,
                coil_pos,
                coil_pos + coilspan / 2,
            ]
            ydata = [
                top_height,
                0,
                -coil_height,
                -coil_height - top_height
            ]
            up_or_down = -arrow_head_length
            arrow_y_pos = coil_height * .88
        else:
            xdata = [
                coil_pos - coilspan / 2,
                coil_pos,
                coil_pos,
                coil_pos - coilspan / 2 + 1,
            ]
            ydata = [
                -coil_height - top_height,
                -coil_height,
                0,
                top_height,
            ]
            up_or_down = arrow_head_length
            arrow_y_pos = coil_height * .12

        ET.SubElement(svg, "rect", {
            "x": f"{coil_pos + 2.5}",
            "y": f"{-coil_height + 1}",
            "width": f"5",
            "height": f"{coil_height - 2}",
            "fill": "lightblue",
        })

        ET.SubElement(svg, "path", {
            "d": f"M {xdata[0]} {ydata[0]} "
                 + " ".join([f"L {x} {y}" for (x, y) in zip(xdata[1:], ydata[1:])]),
            "fill": "none",
            "stroke": f"{coil_color}",
            "stroke-width": ".25px",
            "stroke-linejoin": "round",
            "stroke-linecap": "round",
        })

        arrow_points = [
            (coil_pos, -arrow_y_pos),
            (coil_pos - arrow_head_width / 2, -arrow_y_pos - up_or_down),
            (coil_pos + arrow_head_width / 2, -arrow_y_pos - up_or_down),
        ]

        ET.SubElement(svg, "polygon", {
            "points": " ".join([f"{x},{y}" for (x, y) in arrow_points]),
            "fill": f"{coil_color}",
            "stroke": "none",
        })

        ET.SubElement(svg, "circle", {
            "cx": f"{coil_pos}",
            "cy": f"{-coil_height / 2}",
            "r": ".1em",
            "fill": "white",
        })

        ET.SubElement(svg, "text", {
            "x": f"{coil_pos}",
            "y": f"{-coil_height / 2}",
            "text-anchor": "middle",
            "dominant-baseline": "middle",
            "style": "font-size: .15em; font-family: sans-serif;",
        }).text = str(i + 1)

        if i == data.index(abs(key)) and abs(key) != 2:
            if out:
                x = coil_pos + coilspan / 2 - 1
                name = str(abs(key))
            else:
                x = coil_pos - coilspan / 2 + 1
                name = str(abs(key)) + "'"

            ET.SubElement(svg, "path", {
                "d": f"M {x} {top_height} L {x} {20}",
                "stroke": f"{coil_color}",
                "stroke-width": ".25px",
                "stroke-linecap": "round",
            })

            ET.SubElement(svg, "text", {
                "x": f"{x + 1}",
                "y": f"22",
                "fill": f"{coil_color}",
                "style": "font-size: .25em; font-family: sans-serif;",
            }).text = name

        elif i == len(data) - 1 - data[::-1].index(-abs(key)) and abs(key) != 2:
            if out:
                x = coil_pos + coilspan / 2 - 1
            else:
                x = coil_pos - coilspan / 2 + 1

            ET.SubElement(svg, "path", {
                "d": f"M {x} {top_height} L {x} {20}",
                "stroke": f"{coil_color}",
                "stroke-width": ".25px",
                "stroke-linecap": "round",
            })

            ET.SubElement(svg, "text", {
                "x": f"{x + 1}",
                "y": f"22",
                "fill": f"{coil_color}",
                "style": "font-size: .25em; font-family: sans-serif;",
            }).text = str(abs(key)) + "'"

        elif abs(key) == 2 and (i == data.index(key) or data[data.index(key):].index(key)):
            if out:
                x = coil_pos + coilspan / 2 - 1
                name = str(abs(key))
            else:
                x = coil_pos - coilspan / 2 + 1
                name = str(abs(key)) + "'"

            ET.SubElement(svg, "path", {
                "d": f"M {x} {top_height} L {x} {20}",
                "stroke": f"{coil_color}",
                "stroke-width": ".25px",
                "stroke-linecap": "round",
            })

            ET.SubElement(svg, "text", {
                "x": f"{x + 1}",
                "y": f"22",
                "fill": f"{coil_color}",
                "style": "font-size: .25em; font-family: sans-serif;",
            }).text = name

        else:
            if out:
                xdata = [
                    coil_pos + coilspan / 2 - 1,
                    coil_pos + coilspan / 2 - 1,
                    coil_pos - coilspan / 2,
                ]
            else:
                xdata = [
                    coil_pos - coilspan / 2 + 1,
                    coil_pos - coilspan / 2 + 1,
                    coil_pos + coilspan / 2,
                ]

            ydata = [
                top_height,
                base,
                base
            ]

            ET.SubElement(svg, "path", {
                "d": f"M {xdata[0]} {ydata[0]} " + " ".join([f"L {x} {y}" for (x, y) in zip(xdata[1:], ydata[1:])]),
                "fill": "none",
                "stroke": f"{coil_color}",
                "stroke-width": ".25px",
                "stroke-linejoin": "round",
                "stroke-linecap": "round",
            })

    ET.ElementTree(svg).write(filename)
