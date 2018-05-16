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


def write_fsl_file(machine, basename, inner=False, outer=False):
    model = NewFslRenderer(basename)
    filename = basename + '_' + machine.geom.kind + '.fsl'
    model.render(machine.geom, filename, inner, outer)


def write_main_fsl_file(machine, machine_inner, machine_outer, basename):
    model = NewFslRenderer(basename)
    filename = basename + '.fsl'
    model.render_main(machine, machine_inner, machine_outer, filename)


def symmetry_search(machine,
                    plt,  # plotter
                    kind,
                    symtol=0.0,
                    is_inner=False,
                    is_outer=False,
                    show_plots=True,
                    debug_mode=False,
                    rows=1,
                    cols=1,
                    num=1):
    logger.info("symmetry search for %s", kind)
    machine.clear_cut_lines()
    if show_plots and debug_mode:
        plt.render_elements(machine.geom, Shape,
                            neighbors=True, title=kind)

    if not machine.find_symmetry(symtol):
        logger.info("{}: no symmetry axis found".format(kind))
        plt.add_emptyplot(rows, cols, num, 'no symmetry axis')

        machine_mirror = machine.get_symmetry_mirror()
        machine_slice = machine
    else:
        if show_plots:
            plt.render_elements(machine.geom, Shape,
                                title=kind+' (symmetrylines)',
                                draw_inside=True,
                                rows=rows, cols=cols, num=num, show=False)
        machine_slice = machine.get_symmetry_slice()
        if machine_slice is None:
            machine.kind = kind
            return machine

        machine_mirror = machine_slice.get_symmetry_mirror()

    if machine_mirror is None:
        logger.info("no mirror found")
        machine_ok = machine_slice
    else:
        if show_plots and debug_mode:
            plt.render_elements(machine_mirror.mirror_geom, Shape,
                                title='Mirror of '+kind,
                                rows=rows, cols=cols, num=num, show=True)

        logger.info("mirror found")
        machine_next_mirror = machine_mirror
        while machine_next_mirror is not None:
            machine_mirror = machine_next_mirror
            machine_next_mirror = machine_mirror.get_symmetry_mirror()
            logger.info("another mirror found")
        machine_ok = machine_mirror

    machine_ok.set_minmax_radius()
    # machine_ok.complete_hull(is_inner, is_outer)
    machine_ok.create_auxiliary_lines()
    machine_ok.set_kind(kind)
    return machine_ok


def converter(dxfile,
              rtol=1e-03,
              atol=1e-03,
              mindist=0.01,
              symtol=0.001,
              split=False,
              inner_name='inner',
              outer_name='outer',
              part=(),
              airgap=0.0,
              airgap2=0.0,
              view_only=False,
              show_plots=True,
              write_fsl=False,
              debug_mode=False):
    layers = ()

    basename = os.path.basename(dxfile).split('.')[0]
    logger.info("start reading %s", basename)

    basegeom = Geometry(dxfshapes(dxfile,
                                  mindist=mindist,
                                  layers=layers),
                        rtol=rtol,
                        atol=atol,
                        split=split)
    logger.info("total elements %s", len(basegeom.g.edges()))

    p = PlotRenderer()
    if view_only:
        p.render_elements(basegeom, Shape,
                          neighbors=True,
                          show=True)
        sys.exit(0)

    machine_base = basegeom.get_machine()
    if show_plots:
        p.render_elements(basegeom, Shape,
                          title='Original',
                          with_hull=False,
                          rows=3, cols=2, num=1, show=debug_mode)

    if not machine_base.is_a_machine():
        logger.info("it's Not a Machine!!")
        sys.exit(1)

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
        logger.info("No arc segment found")
        sys.exit(1)

    machine.clear_cut_lines()
    machine.move_to_middle()
    if show_plots and debug_mode:
        p.render_elements(machine.geom, Shape,
                          title='Areas',
                          with_corners=False, show=True)

    if machine.airgap(airgap, airgap2, symtol):
        p.render_elements(machine.geom, Shape,
                          title='Search for airgap failed',
                          with_corners=False, show=True)
        sys.exit(1)

    if show_plots:
        p.render_elements(basegeom, Shape, neighbors=True,
                          title='Original with nodes',
                          rows=3, cols=2, num=2, show=False)

    machine.repair_hull()
    if machine.has_airgap():
        machine_inner = machine.copy(0.0, 2*np.pi, True, True)
        machine_inner = symmetry_search(machine_inner,
                                        p,  # plot
                                        inner_name,
                                        is_inner=True,
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
                                        is_outer=True,
                                        symtol=symtol,
                                        show_plots=show_plots,
                                        rows=3,  # rows
                                        cols=2,  # columns
                                        num=4)   # start num

        machine_inner.sync_with_counterpart(machine_outer)

        machine_inner.search_subregions()
        machine_outer.search_subregions()

        if machine_inner.geom.area_close_to_endangle(2) > 0:
            logger.info("undo mirror of %s", inner_name)
            machine_inner = machine_inner.undo_mirror()
            machine_inner.sync_with_counterpart(machine_outer)
            machine_inner.search_subregions()

        elif machine_outer.geom.area_close_to_endangle(2) > 0:
            logger.info("undo mirror of %s", outer_name)
            machine_outer = machine_outer.undo_mirror()
            machine_inner.sync_with_counterpart(machine_outer)
            machine_outer.search_subregions()

        machine_inner.delete_tiny_elements(mindist)
        machine_outer.delete_tiny_elements(mindist)

        if show_plots:
            p.render_elements(machine_inner.geom, Shape,
                              draw_inside=True, title=inner_name,
                              rows=3, cols=2, num=5, show=False,
                              # with_nodes=True,
                              fill_areas=True)

            p.render_elements(machine_outer.geom, Shape,
                              draw_inside=True, title=outer_name,
                              rows=3, cols=2, num=6, show=False,
                              # with_nodes=True,
                              fill_areas=True)
            p.show_plot()

            # p.render_areas(machine_inner.geom,
            #                with_nodes=True,
            #                single_view=True)
            # p.render_areas(machine_outer.geom)

        if write_fsl:
            write_fsl_file(machine_inner, basename, True, False)
            write_fsl_file(machine_outer, basename, False, True)
            write_main_fsl_file(machine,
                                machine_inner,
                                machine_outer,
                                basename)

    else:
        # No airgap found
        machine = symmetry_search(machine,
                                  p,  # plot
                                  "No_Airgap",
                                  symtol=symtol,
                                  show_plots=show_plots,
                                  rows=3,  # rows
                                  cols=2,  # cols
                                  num=3)   # start num
        if part:
            if part[0] == 'stator':
                machine.geom.search_stator_subregions(part[1])
            else:
                machine.geom.search_rotor_subregions(part[1])

        if show_plots:
            p.render_elements(machine.geom, Shape,
                              draw_inside=True, title='no airgap',
                              rows=3, cols=2, num=5, show=False,
                              fill_areas=True)
            p.show_plot()

        if write_fsl:
            write_fsl_file(machine, basename)

    logger.info("done")
