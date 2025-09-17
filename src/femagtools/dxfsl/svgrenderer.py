# -*- coding: utf-8 -*-
"""
    femagtools.dxfsl.svgrenderer
    ~~~~~~~~~~~~~~~~~~~~~~~~~~--

    renders geometry to svg


"""
import io
import numpy as np
from femagtools.dxfsl.shape import Shape
from femagtools.dxfsl.area import TYPE_TOOTH, TYPE_YOKE
from femagtools.dxfsl.shape import is_Arc, is_Line, is_Circle
from femagtools.dxfsl.functions import alpha_angle
from femagtools import __version__
import logging

logger = logging.getLogger(__name__)


class SvgRenderer(object):
    """renders geom to SVG """

    def __init__(self, name, suffix="", full=False):
        self.basename = name
        self.suffix = suffix
        self.max_radius = 0
        self.head = []
        self.content = []
        self.tail = []
        self.x_list = []
        self.y_list = []
        self.full = full
        self.x_start = 0
        self.y_start = 0
        self.width = 0
        self.height = 0
        self.legend = {}

    def _write_line(self, line, stroke_width=0.03, color='black'):
        self.x_list.append(line.n1[0])
        self.y_list.append(line.n1[1])
        self.x_list.append(line.n2[0])
        self.y_list.append(line.n2[1])

        self.content.append(
            '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="{}" />'.format(
                line.n1[0],
                line.n1[1],
                line.n2[0],
                line.n2[1],
                color,
                stroke_width))

    def _write_arc(self, arc, stroke_width=0.03, color='black'):
        self.x_list.append(arc.n1[0])
        self.y_list.append(arc.n1[1])
        self.x_list.append(arc.n2[0])
        self.y_list.append(arc.n2[1])
        
        self.content.append(
            '<path fill="none" stroke="{}" stroke-width="{}" d="M {} {} A {} {} 0 0 {} {} {}" />'.format(
                color,
                stroke_width,
                arc.n1[0],
                arc.n1[1],
                arc.radius,
                arc.radius,
                arc.turn_lefthand(arc.n1),
                arc.n2[0],
                arc.n2[1]))

    def _write_circle(self, circle, stroke_width=0.03, color='black', dimensions=True):
        if dimensions:
            self.x_list.append(circle.center[0] + circle.radius)
            self.x_list.append(circle.center[0] - circle.radius)
            self.y_list.append(circle.center[1] + circle.radius)
            self.y_list.append(circle.center[1] - circle.radius)
        
        self.content.append(
            '<circle r="{}" cx="{}" cy="{}" fill="none" stroke="{}" stroke-width="{}" />'.format(
                circle.radius,
                circle.center[0],
                circle.center[1],
                color,
                stroke_width))

    def _write_point(self, n, color='blue', radius=0.3):
        self.content.append(
            '<circle r="{}" cx="{}" cy="{}" fill="{}" stroke="none" />'.format(
                radius,
                n[0],
                n[1],
                color))

    def point(self, n, marker, color='blue'):
        self._write_point(n, color=color)

    def _write_element(self, element, stroke_width=0.03, color='black', nodes=False, dimensions=True):
        if is_Arc(element):
            self._write_arc(element, stroke_width=stroke_width, color=color)
            return
        if is_Circle(element):
            self._write_circle(element, stroke_width=stroke_width, color=color, dimensions=dimensions)
            return
        if is_Line(element):
            self._write_line(element, stroke_width=stroke_width, color=color)
            return

    def _write_geometry(self, machine, nodes=False, stroke_width=0.03):
        for e in machine.geom.elements():
            self._write_element(e, stroke_width=stroke_width)

    def _write_symmetry_lines(self, machine):
        for e in machine.geom.cut_lines:
            if e.n1 is None:
                e.n1 = e.p1
                e.n2 = e.p2
            self._write_element(e, color='darkgreen', stroke_width=0.5)

    def _write_airgap_lines(self, machine):
        for e in machine.geom.airgaps:
            if e.n1 is None:
                e.n1 = e.p1
                e.n2 = e.p2
            self._write_element(e, color='red', stroke_width=0.4, dimensions=False)

    def _write_head(self, x_min, x_max, y_min, y_max):
        overhang = 3
        self.x_start = x_min - overhang
        self.width = x_max - x_min + 2 * overhang
        y_length = y_max - y_min
        self.y_start = - (y_max + overhang)
        self.height = y_length + 2 * overhang

        self.head = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<svg viewBox="{} {} {} {}">'.format(self.x_start, self.y_start, self.width, self.height),
            '<g transform="scale(1, -1)">']

    def _write_tail(self, legend=False):
        if self.max_radius > 0.0:
            self.tail.append(  # outer radius
                '<circle r="{}" cx="0" cy="0" fill="none" stroke="black" stroke-width="0.1"/>'.format(self.max_radius))
        self.tail.append(  # center
            '<circle r="0.1" cx="0" cy="0" fill="red" stroke="black" stroke-width="0.05"/>')
        self.tail.append('</g>')
        if legend:
            self._write_legend()
        self.tail.append('</svg>')

    def _write_path(self, area, mirrored=False, stroke_width=0.03, stroke_color='black'):
        self.content += [
            '<path style="fill:{};stroke:{};stroke-width:{}"'.format(
                area.color(),
                stroke_color,
                stroke_width)]

        n0 = None
        for n1, n2, e in area.list_of_elements(tolerant=True):
            self.x_list.append(n2[0])
            self.y_list.append(n2[1])
            if n0 is None:
                n0 = n1
                self.content.append('     d="M {} {}'.format(n1[0], n1[1]))
                
            if is_Arc(e):
                self.content.append('        A {} {} 0 0 {} {} {}'.format(e.radius,
                                                                          e.radius,
                                                                          e.turn_lefthand(n1),
                                                                          n2[0],
                                                                          n2[1]))
            elif is_Line(e):
                self.content.append('        L {} {}'.format(n2[0],
                                                             n2[1]))

        self.content.append('"/>')

    def _write_point_inside(self, geom, area):
        p = area.get_point_inside(geom)
        self._write_point(p, color='magenta', radius=0.2)

    def _write_area(self, machine, area, mirrored=False, part=0, stroke_width=0.03, points=False):
        startangle = machine.startangle
        endangle = machine.endangle
        self._write_path(area, mirrored=mirrored, stroke_width=stroke_width)
        if points:
            self._write_point_inside(machine.geom, area)

        if area.legend():
            self.legend[area.legend()] = area.color()
        mirrored_area = None

        if mirrored:
            axis_m, axis_n = machine.geom.get_axis_m_n(endangle)
            mirrored_area = area.mirror_area(machine.center, axis_m, axis_n, set_nodes=True)
            self._write_path(mirrored_area, stroke_width=stroke_width)
            if points:
                self._write_point_inside(machine.geom, mirrored_area)

        if part == 0:
            return

        if not mirrored:
            alpha = alpha_angle(startangle, endangle)
        else:
            alpha = alpha_angle(startangle, endangle) * 2

        for x in range(1, part):
            area = area.rotate_area(machine.center, alpha, set_nodes=True)
            self._write_path(area, stroke_width=stroke_width)
            if points:
                self._write_point_inside(machine.geom, area)

        if mirrored:
            for x in range(1, part):
                mirrored_area = mirrored_area.rotate_area(machine.center, alpha, set_nodes=True)
                self._write_path(mirrored_area, stroke_width=stroke_width)
                if points:
                    self._write_point_inside(machine.geom, mirrored_area)

    def _write_areas(self, machine, stroke_width=0.03, points=False):
        geom = machine.geom
        self.max_radius = geom.max_radius
        iron = [a for a in geom.area_list if a.is_iron()]
        air = [a for a in geom.area_list if a.is_air()]
        mag = [a for a in geom.area_list if a.is_magnet()]
        wnd = [a for a in geom.area_list if a.is_winding()]

        part = 0
        mirrored = False
        if self.full:
            part = int(machine.get_symmetry_part())
            mirrored = machine.is_mirrored()

        for area in iron:
            self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width, points=points)
        for area in air:
            self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width, points=points)
        for area in mag:
           self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width, points=points)
        for area in wnd:
           self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width, points=points)

    def _write_legend(self):
        if len(self.legend) < 1:
            return

        scalefactor = self.width / 50
        margin = 0.5
        spacing = 0.5
        linespace = 0.3
        colorwidth = 1
        colorheight = 0.7
        fontsize = 1
        chars = max([len(key) for key, value in self.legend.items()])
        strokewidth = 0.1
        ry = 0.5

        margin = margin * scalefactor
        spacing = spacing * scalefactor
        linespace = linespace * scalefactor
        colorwidth = colorwidth * scalefactor
        colorheight = colorheight * scalefactor
        fontsize = fontsize * scalefactor
        strokewidth = strokewidth * scalefactor
        ry = ry * scalefactor

        family = "Helvetica, Arial, sans-serif"
        entries = len(self.legend)

        h = (entries * fontsize) + 2 * spacing + (entries - 1) * linespace
        x = self.x_start + margin
        y = self.y_start + margin
        h = (entries * fontsize) + 2 * spacing + (entries - 1) * linespace
        w = spacing + colorwidth + spacing + fontsize * 0.6 * chars + spacing

        self.tail += [
            '<rect x="{}" y="{}" ry="{}" width="{}" height="{}"'.format(x, y, ry, w, h),
            '      stroke="black" fill="ivory" stroke-width="{}" />'.format(strokewidth)]
        x = x + spacing
        y = y + spacing
        for key, value in self.legend.items():
            self.tail.append(
                '<rect x="{}" y="{}" width="{}" height="{}" stroke="none" fill="{}" />'.format(
                    x,
                    y + (fontsize - colorheight),
                    colorwidth,
                    colorheight,
                    value))
            self.tail.append(
                '<text x="{}" y="{}" font-family="{}" font-size="{}">{}</text>'.format(
                    x + colorwidth + spacing,
                    y + fontsize,
                    family,
                    fontsize,
                    key))
            y = y + fontsize + linespace
        return

    def render(self, machine, no_areas=False, nodes=False, stroke_width=0.03, points=False):
        '''create svg statements with nodechains'''
        if no_areas:
            self._write_geometry(machine, stroke_width=stroke_width)
        else:
            self._write_areas(machine, stroke_width=stroke_width, points=points)

        self._write_symmetry_lines(machine)
        self._write_airgap_lines(machine)
        if nodes:
            machine.geom.render_neighbors(self)

    def write(self, legend=False):
        max_x = max(self.x_list)
        min_x = min(self.x_list)
        max_y = max(self.y_list)
        min_y = min(self.y_list)

        self._write_head(min_x, max_x, min_y, max_y)
        self._write_tail(legend=legend)

        name = self.basename
        if self.suffix:
            name += "-" + self.suffix
        with io.open(name + '.svg', 'w',
                     encoding='utf-8') as f:
            f.write('\n'.join(self.head))
            f.write('\n')
            f.write('\n'.join(self.content))
            f.write('\n')
            f.write('\n'.join(self.tail))
