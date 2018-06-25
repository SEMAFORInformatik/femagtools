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
from femagtools.dxfsl.fslrenderer import FslRenderer
from femagtools.dxfsl.plotrenderer import PlotRenderer
import logging
import logging.config
import numpy as np

logger = logging.getLogger(__name__)


def usage(name):
    print("Usage: ", name,
          " [-h] [--help]")


def write_fsl_file(machine,
                   basename,
                   inner=False,
                   outer=False):
    model = FslRenderer(basename)
    filename = basename + '_' + machine.geom.kind + '.fsl'
    model.render(machine, filename, inner, outer)
    return filename


def write_main_fsl_file(machine,
                        machine_inner,
                        machine_outer,
                        params,
                        basename):
    model = FslRenderer(basename)
    filename = basename + '.fsl'
    model.render_main(machine, machine_inner, machine_outer, params, filename)
    return filename


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
        if show_plots:
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
    conv = {}

    basename = os.path.basename(dxfile).split('.')[0]
    logger.info("start reading %s", basename)

    if part:
        if part[0] not in ('rotor', 'stator'):
            logger.error('FATAL: Parameter rotor or stator expected')
            sys.exit(1)
        if part[1] not in ('in', 'out'):
            logger.error('"{}" has to be defined in/out'.format(part[0]))
            sys.exit(1)

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
            inner_filename = write_fsl_file(machine_inner,
                                            basename,
                                            True,
                                            False)
            outer_filename = write_fsl_file(machine_outer,
                                            basename,
                                            False,
                                            True)

            if machine_inner.geom.is_rotor():
                conv['filename_rotor'] = inner_filename
                conv['filename_stator'] = outer_filename
            else:
                conv['filename_stator'] = inner_filename
                conv['filename_rotor'] = outer_filename

            params = create_femag_parameters(machine_inner,
                                             machine_outer)

            conv['main_filename'] = write_main_fsl_file(machine,
                                                        machine_inner,
                                                        machine_outer,
                                                        params,
                                                        basename)
            conv.update(params)
    else:
        # No airgap found
        name = "No_Airgap"
        inner = False
        outer = False
        params = None

        if part:
            if part[1] == 'in':
                name = inner_name
                inner = True
            else:
                name = outer_name
                outer = True

        machine = symmetry_search(machine,
                                  p,  # plot
                                  name,
                                  symtol=symtol,
                                  show_plots=show_plots,
                                  rows=3,  # rows
                                  cols=2,  # cols
                                  num=3)   # start num
        if part:
            if part[0] == 'stator':
                machine.geom.set_stator()
                machine.geom.search_stator_subregions(part[1])
                params = create_femag_parameters_stator(machine,
                                                        part[1])
            else:
                machine.geom.set_rotor()
                machine.geom.search_rotor_subregions(part[1])
                params = create_femag_parameters_rotor(machine,
                                                       part[1])
        if show_plots:
            p.render_elements(machine.geom, Shape,
                              draw_inside=True, title=name,
                              rows=3, cols=2, num=5, show=False,
                              fill_areas=True)
            p.show_plot()

        if write_fsl:
            conv['filename'] = write_fsl_file(machine, basename, inner, outer)
            if params:
                conv.update(params)

    logger.info("done")
    return conv


def create_femag_parameters(m_inner, m_outer):
    if not (m_inner and m_outer):
        return {}

    params = {}
    geom_inner = m_inner.geom
    geom_outer = m_outer.geom

    parts_inner = int(m_inner.get_symmetry_part())
    parts_outer = int(m_outer.get_symmetry_part())

    if parts_inner > parts_outer:
        num_slots = parts_inner
        num_poles = parts_outer
        num_sl_gen = int(geom_inner.get_symmetry_copies()+1)
        alfa_slot = geom_inner.get_alfa()
        alfa_pole = geom_outer.get_alfa()
    else:
        num_slots = parts_outer
        num_poles = parts_inner
        num_sl_gen = int(geom_outer.get_symmetry_copies()+1)
        alfa_slot = geom_outer.get_alfa()
        alfa_pole = geom_inner.get_alfa()

    params['tot_num_slot'] = num_slots
    params['num_sl_gen'] = num_sl_gen
    params['num_poles'] = num_poles

    params['dy1'] = 2*geom_outer.max_radius
    params['da1'] = 2*geom_outer.min_radius
    params['da2'] = 2*geom_inner.max_radius
    params['dy2'] = 2*geom_inner.min_radius

    params['alfa_slot'] = alfa_slot
    params['alfa_pole'] = alfa_pole
    assert(np.isclose(alfa_slot * num_slots,
                      alfa_pole * num_poles))
    return params


def create_femag_parameters_stator(motor, position):
    params = {}
    num_slots = motor.get_symmetry_part()
    params['tot_num_slot'] = num_slots
    if position == 'in':
        params['da2'] = 2*motor.geom.max_radius
        params['dy2'] = 2*motor.geom.min_radius
    else:
        params['dy1'] = 2*motor.geom.max_radius
        params['da1'] = 2*motor.geom.min_radius
    return params


def create_femag_parameters_rotor(motor, position):
    params = {}
    num_poles = motor.get_symmetry_part()
    params['num_poles'] = num_poles
    if position == 'in':
        params['da2'] = 2*motor.geom.max_radius
        params['dy2'] = 2*motor.geom.min_radius
    else:
        params['dy1'] = 2*motor.geom.max_radius
        params['da1'] = 2*motor.geom.min_radius
    return params
