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

    def _write_line(self, line, stroke_width=0.03, color='black', nodes=False):
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
        if nodes:
            self._write_node(line.n1)
            self._write_node(line.n2)

    def _write_arc(self, arc, stroke_width=0.03, color='black', nodes=False):
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
        if nodes:
            self._write_node(arc.n1)
            self._write_node(arc.n2)

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

    def _write_node(self, n):
        self.content.append(
            '<circle r="0.5" cx="{}" cy="{}" fill="blue" stroke="none" />'.format(
                n[0],
                n[1]))

    def _write_element(self, element, stroke_width=0.03, color='black', nodes=False, dimensions=True):
        if is_Arc(element):
            self._write_arc(element, stroke_width=stroke_width, color=color, nodes=nodes)
            return
        if is_Circle(element):
            self._write_circle(element, stroke_width=stroke_width, color=color, dimensions=dimensions)
            return
        if is_Line(element):
            self._write_line(element, stroke_width=stroke_width, color=color, nodes=nodes)
            return

    def _write_geometry(self, machine, nodes=False, stroke_width=0.03):
        for e in machine.geom.elements():
            self._write_element(e, nodes=nodes, stroke_width=stroke_width)

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
            self._write_element(e, color='red', stroke_width=0.5, dimensions=False)

    def _write_head(self, x_min, x_max, y_min, y_max):
        overhang = 3
        x_start = x_min - overhang
        x_width = x_max - x_min + 2 * overhang
        y_length = y_max - y_min
        y_start = - (y_max + overhang)
        y_height = y_length + 2 * overhang

        self.head = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<svg viewBox="{} {} {} {}">'.format(x_start, y_start, x_width, y_height),
            '<g transform="scale(1, -1)">']

    def _write_tail(self):
        if self.max_radius > 0.0:
            self.tail += [
                '<circle r="{}" cx="0" cy="0" fill="none" stroke="black" stroke-width="0.1"/>'.format(self.max_radius)]
        self.tail += [
            '<circle r="0.1" cx="0" cy="0" fill="red" stroke="black" stroke-width="0.05"/>']
        self.tail += [
            '</g>',
            '</svg>']

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

    def _write_area(self, machine, area, mirrored=False, part=0, stroke_width=0.03):
        startangle = machine.startangle
        endangle = machine.endangle
        self._write_path(area, mirrored=mirrored, stroke_width=stroke_width)
        mirrored_area = None

        if mirrored:
            axis_m, axis_n = machine.geom.get_axis_m_n(endangle)
            mirrored_area = area.mirror_area(machine.center, axis_m, axis_n, set_nodes=True)
            self._write_path(mirrored_area, stroke_width=stroke_width)

        if part == 0:
            return

        if not mirrored:
            alpha = alpha_angle(startangle, endangle)
        else:
            alpha = alpha_angle(startangle, endangle) * 2

        for x in range(1, part):
            area = area.rotate_area(machine.center, alpha, set_nodes=True)
            self._write_path(area, stroke_width=stroke_width)

        if mirrored:
            for x in range(1, part):
                mirrored_area = mirrored_area.rotate_area(machine.center, alpha, set_nodes=True)
                self._write_path(mirrored_area, stroke_width=stroke_width)

    def _write_areas(self, machine, stroke_width=0.03):
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
            self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width)
        for area in air:
            self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width)
        for area in mag:
           self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width)
        for area in wnd:
           self._write_area(machine, area, mirrored=mirrored, part=part, stroke_width=stroke_width)

    def render(self, machine, no_areas=False, nodes=False, stroke_width=0.03):
        '''create svg statements with nodechains'''
        if no_areas:
            self._write_geometry(machine, nodes=nodes, stroke_width=stroke_width)
        else:
            self._write_areas(machine, stroke_width=stroke_width)

        self._write_symmetry_lines(machine)
        self._write_airgap_lines(machine)

    def write(self):
        max_x = max(self.x_list)
        min_x = min(self.x_list)
        max_y = max(self.y_list)
        min_y = min(self.y_list)

        self._write_head(min_x, max_x, min_y, max_y)
        self._write_tail()

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
