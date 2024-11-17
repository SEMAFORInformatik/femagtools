import re
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mpl

rgbpat = re.compile(r'rgb\((\d+),(\d+),(\d+)\)')

def fieldlines(svgfilename, ax=0):
    """plot fieldlines from svg file"""
    lines = []
    cols = []
    tree = ET.parse(svgfilename)
    root = tree.getroot()
    for line in root.findall('svg:line',
                             {'svg': 'http://www.w3.org/2000/svg'}):
        lines.append(
            ((float(line.attrib['x1']), float(line.attrib['x2'])),
             (float(line.attrib['y1']), float(line.attrib['y2']))))
        cols.append(
            [int(x)/256
             for x in rgbpat.findall(line.attrib['stroke'])[0]])
    a = np.array(lines)
    xmax, xmin = np.max(a[:, 0]), np.min(a[:, 0])
    ymax, ymin = np.max(a[:, 1]), np.min(a[:, 1])

    if ax == 0:
        ax = plt.gca()
    ax.set_frame_on(False)
    ax.set_aspect(1)
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymax, ymin))  # upside down
    ax.set_yticks([])
    ax.set_xticks([])

    for l, c in zip(lines, cols):
        ax.add_line(mpl.Line2D(l[0], l[1], color=c))
