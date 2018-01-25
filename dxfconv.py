#!/usr/bin/env python
#
# read a dxf file and create a plot or fsl file
#
# Author: Ronald Tanner
# Date: 2016/01/24
#
import sys
import os
import femagtools.dxfsl.geom as dg
import femagtools.dxfsl.renderer as dr
import argparse
import logging
import logging.config
import numpy as np
# import io

logger = logging.getLogger(__name__)


def usage(name):
    print("Usage: ", name,
          " [-h] [--help]")


def write_fsl(machine, basename, inner=False, outer=False):
    model = dr.NewFslRenderer(basename)
    filename = basename + '_' + machine.geom.kind + '.fsl'
    model.render(machine.geom, filename, inner, outer)


def write_main_fsl(machine, machine_inner, machine_outer, basename):
    model = dr.NewFslRenderer(basename)
    filename = basename + '.fsl'
    model.render_main(machine, machine_inner, machine_outer, filename)


def symmetry_search(machine, kind, sym_tolerance, show_plots,
                    rows=1, cols=1, num=1):
    machine.clear_cut_lines()
    if show_plots and args.debug:
        p.render_elements(machine.geom, dg.Shape, neighbors=True, title=kind)
        # p.render_areas(machine.geom, with_nodes=False, single_view=True)

    if not machine.find_symmetry(sym_tolerance):
        if args.debug:
            print("no symmetry axis found")
        logger.info("{}: no symmetry axis found".format(kind))
        machine_mirror = machine.get_symmetry_mirror()
        machine_slice = machine
    else:
        if show_plots:
            p.render_elements(machine.geom, dg.Shape,
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
        if show_plots and args.debug:
            p.render_elements(machine_mirror.mirror_geom, dg.Shape,
                              title='Mirror of '+kind,
                              rows=rows, cols=cols, num=num, show=True)
        machine_ok = machine_mirror

    machine_ok.complete_hull()
    machine_ok.create_auxiliary_lines()

    if show_plots:
        if args.debug:
            p.render_elements(machine_ok.geom, dg.Shape,
                              draw_inside=True, neighbors=True, title=kind)
            # p.render_areas(machine_ok.geom,
            #                single_view=True, with_nodes=True)
        else:
            p.render_elements(machine_ok.geom, dg.Shape,
                              draw_inside=True, title=kind,
                              rows=rows, cols=cols, num=num+2, show=False)

    machine_ok.set_kind(kind)
    return machine_ok


#############################
#            Main           #
#############################

if __name__ == "__main__":
    loglevel = logging.INFO

    argparser = argparse.ArgumentParser(
        description='Process DXF file and create a plot or FSL file.')
    argparser.add_argument('dxfile',
                           help='name of DXF file')
    argparser.add_argument('-f', '--fsl',
                           help='create fsl',
                           action="store_true")
    argparser.add_argument('--inner',
                           help='name of inner element',
                           dest='inner',
                           default='inner')
    argparser.add_argument('--outer',
                           help='name of outer element',
                           dest='outer',
                           default='outer')
    argparser.add_argument('-a', '--airgap',
                           help='correct airgap',
                           dest='airgap',
                           type=float,
                           default=0.0)
    argparser.add_argument('--airgap2',
                           help='correct airgap',
                           dest='airgap2',
                           type=float,
                           default=0.0)
    argparser.add_argument('-t', '--symtol',
                           help='absolut tolerance to find symmetry axis',
                           dest='sym_tolerance',
                           type=float,
                           default=0.0)
    argparser.add_argument('--rtol',
                           help='relative tolerance (pickdist)',
                           dest='rtol',
                           type=float,
                           default=1e-03)
    argparser.add_argument('--atol',
                           help='absolut tolerance (pickdist)',
                           dest='atol',
                           type=float,
                           default=1e-03)
    argparser.add_argument('-s', '--split',
                           help='split intersections',
                           dest='split',
                           action="store_true")
    argparser.add_argument('-p', '--plot',
                           help='show plots',
                           dest='show_plots',
                           action="store_true")
    argparser.add_argument('-v', '--view',
                           help='show a view only',
                           dest='view',
                           action="store_true")
    argparser.add_argument('-d', '--debug',
                           help='print debug information',
                           dest='debug',
                           action="store_true")
    argparser.add_argument('-l', '--log',
                           help='print debug information',
                           dest='debug',
                           action="store_true")

    args = argparser.parse_args()
    if args.debug:
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel,
                        format='%(asctime)s %(message)s')
    if args.airgap > 0.0:
        if args.debug:
            if args.airgap2 > 0.0:
                print("Airgap is set from {} to {}".
                      format(args.airgap, args.airgap2))
            else:
                print("Airgap is set to {}".format(args.airgap))
        else:
            if args.airgap2 > 0.0:
                logger.info("Airgap is set from {} to {}".
                            format(args.airgap, args.airgap2))
            else:
                logger.info("Airgap is set to {}".format(args.airgap))

    layers = ()

    incl_bnd = True

    rtol = args.rtol
    atol = args.atol

    inner_name = args.inner
    outer_name = args.outer

    basename = os.path.basename(args.dxfile).split('.')[0]
    logger.info("start reading %s", basename)

    basegeom = dg.Geometry(dg.dxfshapes(args.dxfile, layers=layers),
                           rtol=rtol,
                           atol=atol,
                           split=args.split)
    logger.info("total elements %s", len(basegeom.g.edges()))

    p = dr.PlotRenderer()
    if args.view:
        p.render_elements(basegeom, dg.Shape, neighbors=True)
        sys.exit(0)

    machine_base = basegeom.get_machine()
    if args.show_plots:
        p.render_elements(basegeom, dg.Shape,
                          title='Original',
                          rows=3, cols=2, num=1, show=args.debug)

    if not machine_base.is_a_machine():
        print("it's Not a Machine!!")
        sys.exit(1)

    if not machine_base.is_full():
        machine_base.repair_hull()
        if args.show_plots and args.debug:
            print("===== Original (REPAIRED HULL) =====")
            p.render_elements(basegeom, dg.Shape, with_corners=True, show=True)

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
    if args.show_plots and args.debug:
        print("===== Areas =====")
        # p.render_areas(machine.geom, with_nodes=True)
        p.render_elements(machine.geom, dg.Shape,
                          with_corners=False, show=True)

    machine.airgap(args.airgap, args.airgap2, args.sym_tolerance)

    if machine.has_airgap():
        if args.show_plots:
            p.render_elements(basegeom, dg.Shape, neighbors=True,
                              title='Original with nodes',
                              rows=3, cols=2, num=2, show=False)

        machine_inner = machine.copy(0.0, 2*np.pi, True, True)
        machine_inner = symmetry_search(machine_inner, inner_name,
                                        args.sym_tolerance,
                                        args.show_plots, 3, 2, 3)
        machine_inner.set_inner()

        machine_outer = machine.copy(0.0, 2*np.pi, True, False)
        machine_outer = symmetry_search(machine_outer, outer_name,
                                        args.sym_tolerance,
                                        args.show_plots, 3, 2, 4)

        machine_inner.sync_with_counterpart(machine_outer)
        p.show_plot()

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

        if args.fsl:
            write_fsl(machine_inner, basename, True, False)
            write_fsl(machine_outer, basename, False, True)
            write_main_fsl(machine, machine_inner, machine_outer, basename)

    else:
        machine = symmetry_search(machine, "No_Airgap",
                                  args.sym_tolerance, args.show_plots, 4, 2, 3)
        p.show_plot()

        if args.fsl:
            machine.search_subregions()
            write_fsl(machine, basename)

    logger.info("done")
