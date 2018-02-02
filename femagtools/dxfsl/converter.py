#!/usr/bin/env python
#
# read a dxf file and create a plot or fsl file
#
# Author: Ronald Tanner
# Date: 2016/01/24
#
import sys
import os
from femagtools.dxfsl.geom import Geometry, dxfshapes
from femagtools.dxfsl.shape import Shape
from femagtools.dxfsl.renderer import NewFslRenderer, PlotRenderer
import logging
import logging.config
import numpy as np

logger = logging.getLogger(__name__)


def usage(name):
    print("Usage: ", name,
          " [-h] [--help]")


def write_fsl(machine, basename, inner=False, outer=False):
    model = NewFslRenderer(basename)
    filename = basename + '_' + machine.geom.kind + '.fsl'
    model.render(machine.geom, filename, inner, outer)


def write_main_fsl(machine, machine_inner, machine_outer, basename):
    model = NewFslRenderer(basename)
    filename = basename + '.fsl'
    model.render_main(machine, machine_inner, machine_outer, filename)


def symmetry_search(machine,
                    plt,  # plotter
                    kind,
                    symtol=0.0,
                    show_plots=True,
                    write_fsl=False,
                    debug_mode=False,
                    rows=1,
                    cols=1,
                    num=1):
    machine.clear_cut_lines()
    if show_plots and debug_mode:
        plt.render_elements(machine.geom, Shape,
                            neighbors=True, title=kind)

    if not machine.find_symmetry(symtol):
        if debug_mode:
            print("no symmetry axis found")
        logger.info("{}: no symmetry axis found".format(kind))
        machine_mirror = machine.get_symmetry_mirror()
        machine_slice = machine
    else:
        if show_plots:
            plt.render_elements(machine.geom, Shape,
                                title=kind+' (symmetrylines)',
                                rows=rows, cols=cols, num=num, show=False)
        machine_slice = machine.get_symmetry_slice()
        if machine_slice is None:
            machine.kind = kind
            return machine

        machine_mirror = machine_slice.get_symmetry_mirror()

    if machine_mirror is None:
        machine_ok = machine_slice
    else:
        if show_plots and debug_mode:
            plt.render_elements(machine_mirror.mirror_geom, Shape,
                                title='Mirror of '+kind,
                                rows=rows, cols=cols, num=num, show=True)
        machine_ok = machine_mirror

    machine_ok.complete_hull()
    machine_ok.create_auxiliary_lines()

    machine_ok.set_kind(kind)
    return machine_ok


def converter(dxfile,
              rtol=1e-03,
              atol=1e-03,
              symtol=0.0,
              split=False,
              inner_name='inner',
              outer_name='outer',
              airgap=0.0,
              airgap2=0.0,
              view_only=False,
              show_plots=True,
              debug_mode=False):
    layers = ()

    basename = os.path.basename(dxfile).split('.')[0]
    logger.info("start reading %s", basename)

    basegeom = Geometry(dxfshapes(dxfile, layers=layers),
                        rtol=rtol,
                        atol=atol,
                        split=split)
    logger.info("total elements %s", len(basegeom.g.edges()))

    p = PlotRenderer()
    if view_only:
        p.render_elements(basegeom, Shape, neighbors=True)
        sys.exit(0)

    machine_base = basegeom.get_machine()
    if show_plots:
        p.render_elements(basegeom, Shape,
                          title='Original',
                          rows=3, cols=2, num=1, show=debug_mode)

    if not machine_base.is_a_machine():
        print("it's Not a Machine!!")
        sys.exit(1)

    if not machine_base.is_full():
        machine_base.repair_hull()
        if show_plots and debug_mode:
            print("===== Original (REPAIRED HULL) =====")
            p.render_elements(basegeom, Shape, with_corners=True, show=True)

    if machine_base.is_full() or \
       machine_base.is_half() or \
       machine_base.is_quarter():
        # create a copy for further processing
        machine = machine_base.full_copy()
    else:
        # machine shape is unclear
        machine_base.set_center(0.0, 0.0)
        machine_base.set_radius(9999999)
        machine = machine_base.full_copy()

    if machine.part_of_circle() == 0:
        print("No arc segment found")
        sys.exit(1)

    machine.clear_cut_lines()
    machine.move_to_middle()
    if show_plots and debug_mode:
        print("===== Areas =====")
        p.render_elements(machine.geom, Shape,
                          with_corners=False, show=True)

    machine.airgap(airgap, airgap2, symtol)

    if machine.has_airgap():
        if show_plots:
            p.render_elements(basegeom, Shape, neighbors=True,
                              title='Original with nodes',
                              rows=3, cols=2, num=2, show=False)

        machine_inner = machine.copy(0.0, 2*np.pi, True, True)
        machine_inner = symmetry_search(machine_inner,
                                        p,  # plot
                                        inner_name,
                                        symtol=symtol,
                                        show_plots=show_plots,
                                        rows=3,  # rows
                                        cols=2,  # columns
                                        num=3)   # start num
        machine_inner.set_inner()

        machine_outer = machine.copy(0.0, 2*np.pi, True, False)
        machine_outer = symmetry_search(machine_outer,
                                        p,  # plot
                                        outer_name,
                                        symtol=symtol,
                                        show_plots=show_plots,
                                        rows=3,  # rows
                                        cols=2,  # columns
                                        num=4)   # start num

        machine_inner.sync_with_counterpart(machine_outer)

        machine_inner.search_subregions()
        machine_outer.search_subregions()

        if machine_inner.geom.area_close_to_endangle(2) > 0:
            machine_inner.undo_mirror()
            machine_inner.sync_with_counterpart(machine_outer)
            machine_inner.search_subregions()

        elif machine_outer.geom.area_close_to_endangle(2) > 0:
            machine_outer.undo_mirror()
            machine_inner.sync_with_counterpart(machine_outer)
            machine_outer.search_subregions()

        p.render_elements(machine_inner.geom, Shape,
                          draw_inside=True, title=inner_name,
                          rows=3, cols=2, num=5, show=False,
                          fill_areas=False)

        p.render_elements(machine_outer.geom, Shape,
                          draw_inside=True, title=outer_name,
                          rows=3, cols=2, num=6, show=False,
                          fill_areas=False)
        p.show_plot()

        if write_fsl:
            write_fsl(machine_inner, basename, True, False)
            write_fsl(machine_outer, basename, False, True)
            write_main_fsl(machine, machine_inner, machine_outer, basename)

    else:
        machine = symmetry_search(machine,
                                  p,  # plot
                                  "No_Airgap",
                                  symtol=symtol,
                                  show_plots=show_plots,
                                  rows=3,  # rows
                                  cols=2,  # cols
                                  num=3)   # start num
        p.show_plot()

        if write_fsl:
            machine.search_subregions()
            write_fsl(machine, basename)

    logger.info("done")
