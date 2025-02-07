"""read a dxf file and create a plot or fsl file

"""
import os
import io
from pathlib import Path
from femagtools import __version__
from femagtools.dxfsl.geom import Geometry
from femagtools.dxfsl.shape import Shape
from femagtools.dxfsl.fslrenderer import FslRenderer, agndst
from femagtools.dxfsl.plotrenderer import PlotRenderer
from femagtools.dxfsl.concat import Concatenation
from femagtools.dxfsl.functions import Timer, middle_angle
from femagtools.dxfsl.journal import Journal, getJournal
from femagtools.dxfsl.area import TYPE_WINDINGS
from femagtools.dxfsl.areabuilder import disable_logging, enable_logging
import logging
import logging.config
import numpy as np
import sys

logger = logging.getLogger(__name__)


def plot_geom(doit, plt, geom, title="Plot", areas=True):
    if not doit or not plt:
        return

    logger.info("Prepare Plot %s", title)
    plt.render_elements(geom, Shape,
                        draw_inside=areas,
                        draw_groups=False,
                        draw_phi=True,
                        critical=True,
                        title=title,
                        show=True,
                        with_corners=True,
                        with_nodes=False,
                        neighbors=True,
                        write_id=areas,
                        with_legend=False,
                        fill_areas=areas)


def symmetry_search(machine,
                    plt=None,  # plotter
                    kind="single",
                    mindist=0.01,
                    symtol=0.0,
                    sympart=0,
                    is_inner=False,
                    is_outer=False,
                    show_plots=True,
                    debug_mode=False,
                    rows=1,
                    cols=1,
                    num=1):
    logger.info("*** Begin symmetry search for %s ***", kind)

    def return_machine(machine, kind):
        machine.set_attributes(kind=kind, inner=is_inner, outer=is_outer)
        machine.check_and_correct_geom(kind)
        machine.delete_tiny_elements(mindist)
        return machine

    machine.set_attributes(kind=kind, inner=is_inner, outer=is_outer)
    machine.clear_cut_lines()
    if show_plots and debug_mode:
        plt.render_elements(machine.geom, Shape,
                            neighbors=True, title=kind)

    if sympart > 0:
        if not machine.is_full():
            logger.error("force symmetry failed")
            sys.exit(1)
        machine_ok = machine.get_forced_symmetry(sympart)
        logger.info("*** End of symmetry search for %s ***", kind)
        return return_machine(machine_ok, kind)

    plot_geom(False,  # for developer
              plt, machine.geom,
              title="Before Symmetry ({})".format(kind))

    if machine.find_symmetry(symtol, is_inner, is_outer, None):
        logger.info(" - {}: symmetry axis found".format(kind))
        plot_geom(False,  # for developer
                  plt, machine.geom,
                  title="Symmetry found")

        if show_plots:
            plt.render_elements(machine.geom, Shape,
                                title=kind+' (symmetrylines)',
                                draw_inside=True,
                                rows=rows, cols=cols, num=num, show=False)
        machine_slice = machine.get_symmetry_slice()
        if machine_slice is None:
            logger.info(" - no slice extracted ?!?")
            return return_machine(machine, kind)

        plot_geom(False,  # for developer
                  plt, machine_slice.geom,
                  title="Slice of {}".format(kind))

        # no third after slice
        machine_mirror = machine_slice.get_symmetry_mirror(no_third=True)

        if machine_mirror is not None:
            logger.info("*** End of symmetry search for %s (symmetry and mirror) ***", kind)
            return return_machine(machine_mirror, kind)

        logger.info("*** End of symmetry search for %s (symmetry) ***", kind)
        return return_machine(machine_slice, kind)

    # --- no symmetry slice found ---
    logger.info(" - {}: no symmetry axis found".format(kind))
    if show_plots:
        if debug_mode:
            plt.render_elements(machine.geom, Shape,
                                title=kind+' (no symmetry axis)',
                                draw_inside=True,
                                neighbors=True,
                                rows=rows, cols=cols, num=num, show=False)
        else:
            plt.add_emptyplot(rows, cols, num, 'no symmetry axis')

    if not machine.is_startangle_zero():
        logger.debug("Rotate geometry to 0.0")
        machine.rotate_to(0.0)
        machine.set_alfa_and_corners()

    plot_geom(False,  # for developer
              plt, machine.geom,
              title="Before Mirror ({})".format(kind))

    machine_mirror = machine.get_symmetry_mirror()
    machine_slice = machine
    machine_slice.set_alfa_and_corners()

    if machine_mirror is None:
        logger.info(" - no mirror found")
        machine_ok = machine_slice
    else:
        plot_geom(False,  # for developer
                  plt, machine_mirror.geom,
                  title="After Mirror ({})".format(kind))

        logger.info(" - mirror found")
        machine_next_mirror = machine_mirror.get_symmetry_mirror()
        while machine_next_mirror is not None:
            logger.info(" - another mirror found")
            machine_mirror = machine_next_mirror
            plot_geom(False,  # for developer
                      plt, machine_mirror.geom,
                      title="Next Mirror ({})".format(kind))
            machine_next_mirror = machine_mirror.get_symmetry_mirror()

        machine_ok = machine_mirror

    if not machine_ok.is_startangle_zero():
        logger.debug("Rotate geometry to 0.0")
        machine_ok.rotate_to(0.0)
        machine_ok.set_alfa_and_corners()

    logger.info("*** End of symmetry search for %s ***", kind)
    return return_machine(machine_ok, kind)


def build_machine_rotor(machine, inner, mindist, plt, EESM=False, single=False):
    logger.debug("Begin of build_machine_rotor")

    if machine.has_windings():
        logger.debug("do nothing here with windings in rotor")
        logger.debug("End of build_machine_rotor")
        return machine

    timer = Timer(start_it=True)

    if machine.is_mirrored():
        logger.debug("Rotor is mirrored")
        machine_temp = machine.undo_mirror()
        machine_temp.delete_tiny_elements(mindist)
        machine_temp.geom.set_rotor()
        machine_temp.search_rotor_subregions(EESM, single=single)
    else:
        machine_temp = machine

    plot_geom(False,  # for developer
              plt, machine_temp.geom,
              title="Inner Rotor check magnets")

    machine_slice = machine_temp.get_forced_magnet_slice()
    if machine_slice:
        plot_geom(False,  # for developer
                  plt, machine_slice.geom,
                  title="Rotor Magnet Slice")

        machine_temp = machine_slice
        machine_temp.geom.set_rotor()
        machine_temp.rebuild_subregions(EESM, single=single)
        plot_geom(False,  # for developer
                  plt, machine_temp.geom,
                  title="Rotor Magnet Slice after Rebuild")

    rebuild = False
    if machine_temp.has_magnets_in_the_middle():
        logger.debug("Magnets cut")
        rebuild = machine_temp.create_mirror_lines_outside_magnets()
    else:
        if machine.is_mirrored():
            logger.debug("Back to the mirrored machine")
            machine_temp = machine  # undo
            rebuild = machine_temp.create_inner_auxiliary_arcs()
        else:
            rebuild = machine_temp.create_mirror_lines_outside_magnets()
    if rebuild:
        machine_temp.geom.create_list_of_areas(delete=True)

    if machine_temp.create_auxiliary_lines():
        logger.debug("Auxiliary Lines created: rebuild subregions")
        rebuild = True
    if rebuild:
        machine_temp.rebuild_subregions(EESM, single=single)

    machine_temp.geom.recalculate_magnet_orientation()
    if inner:
        machine_temp.create_inner_corner_areas()

    if not machine_temp.is_mirrored():
        plot_geom(False,  # for developer
                  plt, machine_temp.geom,
                  title="Rotor before Boundery Corr")
        machine_temp.create_boundary_nodes()

    plot_geom(False,  # for developer
              plt, machine_temp.geom,
              title="Final Rotor")

    t = timer.stop("-- rotor created in %0.4f seconds --")
    journal.put('time_rotor_created', t)

    logger.debug("End of build_machine_rotor")
    return machine_temp


def build_machine_stator(machine, inner, mindist, plt, EESM=False, single=False):
    logger.debug("Begin of build_machine_stator")
    timer = Timer(start_it=True)

    if not machine.geom.is_stator():
        logger.debug("Rotor with windings")

    if machine.is_mirrored():
        plot_geom(False,  # for developer
                  plt, machine.previous_machine.geom,
                  title="Mirrored Stator")

        logger.debug("undo mirrored windings")
        machine_temp = machine.undo_mirror()
        machine_temp.delete_tiny_elements(mindist)
        machine_temp.geom.set_stator()
        machine_temp.search_stator_subregions(single=single)
        if not machine_temp.has_windings_in_the_middle():
            logger.debug("Back to the mirrored machine")
            machine_temp = machine  # undo
        else:
            machine_temp.create_mirror_lines_outside_windings()
    else:
        machine_temp = machine

    if machine_temp.geom.reduce_element_nodes(mindist):
        machine_temp.rebuild_subregions(EESM, single=single)
        plot_geom(False,  # for developer
                  plt, machine_temp.geom,
                  title="Nodes reduced")

    machine_slice = machine_temp.get_forced_winding_slice()
    if machine_slice:
        plot_geom(False,  # for developer
                  plt, machine_slice.geom,
                  title="Stator Winding Slice")

        machine_temp = machine_slice
        machine_temp.geom.set_stator()
        machine_temp.rebuild_subregions(EESM, single=single)
        if machine_temp.has_windings_in_the_middle():
            machine_temp.create_mirror_lines_outside_windings()

        plot_geom(False,  # for developer
                  plt, machine_temp.geom,
                  title="Stator Winding Slice after Rebuild")

    if machine_temp.create_auxiliary_lines():
        machine_temp.rebuild_subregions(EESM, single=single)
        plot_geom(False,  # for developer
                  plt, machine_temp.geom,
                  title="Stator with Auxiliary Lines")

    if inner:
        machine_temp.create_inner_corner_areas()

    if not machine_temp.is_mirrored():
        plot_geom(False,  # for developer
                  plt, machine_temp.geom,
                  title="Stator before Boundary Corr")
        machine_temp.create_boundary_nodes()

    plot_geom(False,  # for developer
              plt, machine_temp.geom,
              title="Final Stator")

    t = timer.stop("-- stator created in %0.4f seconds --")
    journal.put('time_stator_created', t)

    logger.debug("End of build_machine_stator")
    return machine_temp


def build_inner_machine(machine,
                        mindist=0.01,
                        plt=None,
                        EESM=False):
    logger.info("Begin of build_inner_machine")
    machine.search_subregions(EESM)

    if machine.geom.is_rotor():  # Inner mirrored rotor
        machine = build_machine_rotor(machine,
                                      True,  # is inner
                                      mindist,
                                      plt,
                                      EESM=EESM)

    if machine.geom.is_stator() or machine.has_windings():
        machine = build_machine_stator(machine,
                                       True,
                                       mindist,
                                       plt,
                                       EESM=EESM)

    machine.search_critical_elements(mindist)
    logger.info("End of build_inner_machine")
    return machine


def build_outer_machine(machine,
                        mindist=0.01,
                        plt=None,
                        EESM=False):
    logger.info("Begin of build_outer_machine")
    machine.search_subregions(EESM)

    if machine.geom.is_rotor():  # Outer mirrored rotor
        machine = build_machine_rotor(machine,
                                      False,  # is outer
                                      mindist,
                                      plt,
                                      EESM=EESM)

    if machine.geom.is_stator() or machine.has_windings():
        machine = build_machine_stator(machine,
                                       False,
                                       mindist,
                                       plt,
                                       EESM=EESM)

    machine.search_critical_elements(mindist)
    logger.info("End of build_outer_machine")
    return machine


def convert(dxfile,
            EESM=False,
            rtol=1e-04,
            atol=1e-03,
            mindist=0.01,
            symtol=0.001,
            sympart=0,
            split=False,
            inner_name='inner',
            outer_name='outer',
            part=(),
            airgap=0.0,
            airgap2=0.0,
            da=0.0,
            dy=0.0,
            nodedist=1,
            view_only=False,
            view_korr=False,
            show_plots=False,
            show_areas=False,
            small_plots=False,
            write_fsl=True,
            write_fsl_single=False,
            write_png=False,
            write_id=False,
            full_model=False,
            debug_mode=False,
            write_journal=False):
    global journal
    layers = ()
    conv = {}

    input_file = Path(dxfile)
    if not input_file.is_file():
        logger.error("File %s is not available", input_file)
        sys.exit(1)

    basename = input_file.stem
    if part:
        logger.info("***** start processing %s (%s) [%s] *****",
                    basename,
                    part,
                    __version__)
    else:
        logger.info("***** start processing %s [%s] *****",
                    basename,
                    __version__)
    timer = Timer(start_it=True)
    start_timer = Timer(start_it=True)

    journal = getJournal(name='converter_journal', aktiv=write_journal)
    journal.get_journal(input_file.name)
    journal.set_filename(str(input_file.resolve()))
    journal.set('success', False)
    journal.write_journal()

    if part:
        if part[0] not in ('rotor', 'stator'):
            logger.error('FATAL: Parameter rotor or stator expected')
            return dict(error='unknown part {}'.format(part))
        if part[1] not in ('in', 'out'):
            logger.error('"{}" has to be defined in/out'.format(part[0]))
            return dict(error='unknown location {}'.format(part[1]))
    else:
        if da:
            logger.warn("distance airgap (da) ignored")
            da = 0.0
        if dy:
            logger.warn("distance yoke (dy) ignored")
            dy = 0.0

    split_ini = split
    split_cpy = False
    if not (part or view_only):
        split_ini = False
        split_cpy = split

    if write_fsl_single:
        write_fsl = True
    if small_plots:
        show_plots = False

    try:
        if input_file.suffix in ['.fem', '.FEM']:
            from .femparser import femshapes
            basegeom = Geometry(femshapes(dxfile),
                                rtol=rtol,
                                atol=atol,
                                split=split_ini,
                                concatenate=True,
                                connect=True,
                                delete=True,
                                adjust=True,
                                main=True)
        elif input_file.suffix in ['.dxf', '.DXF']:
            from .dxfparser import dxfshapes
            basegeom = Geometry(dxfshapes(dxfile,
                                          mindist=mindist,
                                          layers=layers),
                                rtol=rtol,
                                atol=atol,
                                split=split_ini,
                                concatenate=True,
                                connect=True,
                                delete=True,
                                adjust=True,
                                main=True)
        elif input_file.suffix in ['.svg', '.SVG']:
            from .svgparser import svgshapes
            basegeom = Geometry(svgshapes(dxfile),
                                rtol=rtol,
                                atol=atol,
                                split=split_ini,
                                concatenate=True,
                                connect=True,
                                delete=True,
                                adjust=True,
                                main=True)
        else:
            logger.error("Unexpected file %s", input_file)
            sys.exit(1)
    except FileNotFoundError as ex:
        logger.error(ex)
        return dict()

    logger.info("total elements %s", len(basegeom.g.edges()))

    p = PlotRenderer()

    if view_only:
        logger.info("View only")
        p.render_elements(basegeom, Shape,
                          neighbors=True,
                          png=write_png,
                          show=True)
        return dict()

    plot_geom(False,  # for developer
              p, basegeom,
              title="Before finding Machine")

    machine_base = basegeom.get_machine()
    if show_plots:
        p.render_elements(basegeom, Shape,
                          title=input_file.name,
                          with_hull=False,
                          rows=3, cols=2, num=1, show=debug_mode)

    if not machine_base.is_a_machine():
        logger.warn("it's Not a Machine")
        return dict(error='machine not detected')

    if not (machine_base.part > 0):
        # machine shape is unclear
        logger.warn("machine shape is unclear")
        machine_base.set_center(0.0, 0.0)
        machine_base.set_radius(9999999)

    machine = machine_base

    if machine.part_of_circle() == 0:
        logger.warn("No arc segment found")
        return dict(error='no arc segment found')

    machine.clear_cut_lines()
    machine.move_to_middle()
    if show_plots and debug_mode:
        p.render_elements(machine.geom, Shape,
                          title='Areas',
                          neighbors=True,
                          with_corners=False, show=True)

    if machine.airgap(airgap, airgap2, symtol):
        p.render_elements(machine.geom, Shape,
                          title='Search for airgap failed',
                          with_corners=False, show=True)
        return dict(error='no airgap found')

    if show_plots:
        p.render_elements(basegeom, Shape, neighbors=True,
                          title='Original with nodes',
                          rows=3, cols=2, num=2, show=False)

    machine.repair_hull()

    if machine.has_airgap():
        logger.info("=== airgap is %s ===", machine.airgap_radius)
        machine_inner = machine.copy(startangle=0.0,
                                     endangle=2*np.pi,
                                     airgap=True,
                                     inside=True,
                                     split=split_cpy,
                                     delete_appendices=True,
                                     concatenate=True,
                                     connect=True)
        machine_inner.set_inner()
        machine_outer = machine.copy(startangle=0.0,
                                     endangle=2*np.pi,
                                     airgap=True,
                                     inside=False,
                                     split=split_cpy,
                                     delete_appendices=True,
                                     concatenate=True,
                                     connect=True)
        machine_outer.set_outer()

        # check airgap
        if machine.has_alternative_airgap():
            bad_inner = machine_inner.check_airgap()
            bad_outer = False
            if not bad_inner:
                bad_outer = machine_outer.check_airgap()
            if bad_inner or bad_outer:
                logger.warning("***** bad Airgap *****")
                if bad_inner:
                    plot_geom(False,  # for developer
                              p, machine_inner.geom,
                              title="Bad Inner Airgap Geometry")
                if bad_outer:
                    plot_geom(False,  # for developer
                              p, machine_outer.geom,
                              title="Bad Outer Airgap Geometry")

                if machine.install_alternative_airgap():
                    logger.info("=== alternative airgap is %s ===",
                                machine.airgap_radius)
                    machine_inner = machine.copy(startangle=0.0,
                                                 endangle=2*np.pi,
                                                 airgap=True,
                                                 inside=True,
                                                 split=split_cpy,
                                                 delete_appendices=True,
                                                 concatenate=True,
                                                 connect=True)
                    machine_inner.set_inner()
                    machine_outer = machine.copy(startangle=0.0,
                                                 endangle=2*np.pi,
                                                 airgap=True,
                                                 inside=False,
                                                 split=split_cpy,
                                                 delete_appendices=True,
                                                 concatenate=True,
                                                 connect=True)
                    machine_outer.set_outer()

        start_timer.stop("-- first part in %0.4f seconds --", info=True)

        process_timer = Timer(start_it=True)
        # inner part
        machine_inner = symmetry_search(machine=machine_inner,
                                        plt=p,  # plot
                                        kind=inner_name,
                                        is_inner=True,
                                        mindist=mindist,
                                        symtol=symtol,
                                        show_plots=show_plots,
                                        rows=3,  # rows
                                        cols=2,  # columns
                                        num=3)   # start num

        # outer part
        machine_outer = symmetry_search(machine_outer,
                                        plt=p,  # plot
                                        kind=outer_name,
                                        is_outer=True,
                                        mindist=mindist,
                                        symtol=symtol,
                                        show_plots=show_plots,
                                        rows=3,  # rows
                                        cols=2,  # columns
                                        num=4)   # start num

        process_timer.stop("-- symmetry search in %0.4f seconds --", info=True)

        machine_inner.sync_with_counterpart(machine_outer)

        final_timer = Timer(start_it=True)
        machine_inner = build_inner_machine(machine_inner,
                                            mindist=mindist,
                                            plt=p,
                                            EESM=EESM)

        machine_outer = build_outer_machine(machine_outer,
                                            mindist,
                                            p,
                                            EESM=EESM)
        final_timer.stop("-- final part in %0.4f seconds --", info=True)

        machine_inner.sync_with_counterpart(machine_outer)
        logger.info("***** END of work: %s *****", basename)

        if machine_inner.geom.is_rotor():
            inner_title = "Rotor"
            outer_title = "Stator"
        else:
            inner_title = "Stator"
            outer_title = "Rotor"

        plot_geom(False,  # for developer
                  p, machine_inner.geom,
                  title="Final Inner Geometry")

        plot_geom(False,  # for developer
                  p, machine_outer.geom,
                  title="Final Outer Geometry")

        if show_plots:
            p.render_elements(machine_inner.geom, Shape,
                              draw_inside=True, title=inner_title,
                              rows=3, cols=2, num=5, show=False,
                              with_corners=False,
                              with_nodes=False,
                              neighbors=False,
                              write_id=write_id,
                              draw_phi=True,
                              fill_areas=True)

            p.render_elements(machine_outer.geom, Shape,
                              draw_inside=True, title=outer_title,
                              rows=3, cols=2, num=6, show=False,
                              with_corners=False,
                              with_nodes=False,
                              neighbors=False,
                              write_id=write_id,
                              draw_phi=True,
                              fill_areas=True)
        elif small_plots:
            #p.figure(figsize=(9, 5)).suptitle(input_file.name, fontsize=16)
            p.render_elements(machine_inner.geom, Shape,
                              draw_inside=True, title=inner_title,
                              rows=1, cols=2, num=1, show=False,
                              with_corners=False,
                              with_nodes=False,
                              neighbors=False,
                              write_id=write_id,
                              draw_phi=True,
                              fill_areas=True)

            p.render_elements(machine_outer.geom, Shape,
                              draw_inside=True, title=outer_title,
                              rows=1, cols=2, num=2, show=False,
                              with_corners=False,
                              with_nodes=False,
                              neighbors=False,
                              write_id=write_id,
                              draw_phi=True,
                              fill_areas=True)

        if show_plots or small_plots:
            if write_png:
                p.write_plot(basename)
            else:
                p.show_plot()

        if show_areas:
            p.render_elements(machine_inner.geom, Shape,
                              title=inner_name,
                              show=True,
                              draw_inside=True,
                              neighbors=True,
                              fill_areas=True)
            p.render_areas(machine_inner.geom,
                           title=inner_name,
                           with_nodes=True,
                           single_view=True)

            p.render_elements(machine_outer.geom, Shape,
                              title=outer_name,
                              show=True,
                              draw_inside=True,
                              neighbors=True,
                              fill_areas=True)
            p.render_areas(machine_outer.geom,
                           title=outer_name,
                           with_nodes=True,
                           single_view=True)

        params = create_femag_parameters(machine_inner,
                                         machine_outer,
                                         nodedist)

        if write_fsl:
            if machine_inner.is_full() or machine_outer.is_full():
                logger.warning("it's not possible to create fsl-file")
                return None

            mtype = 'EESM' if EESM else 'PMSM'
            fslrenderer = FslRenderer(basename, mtype)
            inner = fslrenderer.render(machine_inner, inner=True)
            outer = fslrenderer.render(machine_outer, outer=True)
            if write_fsl_single:
                inner_single = fslrenderer.render(machine_inner, inner=True,
                                                  standalone=True)
                outer_single = fslrenderer.render(machine_outer, outer=True,
                                                  standalone=True)
            else:
                inner_single = None
                outer_single = None
            if full_model:
                params['num_sl_gen'] = params.get('tot_num_slot', 0)
            params['agndst'] = agndst(params.get('da1'),
                                      params.get('da2'),
                                      params.get('tot_num_slot'),
                                      params.get('num_poles'),
                                      params.get('nodedist'))

            if params['external_rotor']:
                conv['fsl_rotor'] = outer
                if outer_single:
                    conv['fsl_rotor_single'] = outer_single
                conv['fsl_stator'] = inner
                if inner_single:
                    conv['fsl_stator_single'] = inner_single
            else:
                conv['fsl_rotor'] = inner
                if inner_single:
                    conv['fsl_rotor_single'] = inner_single
                conv['fsl_stator'] = outer
                if outer_single:
                    conv['fsl_stator_single'] = outer_single

            conv['fsl'] = fslrenderer.render_main(
                machine_inner, machine_outer,
                inner, outer,
                params)
    else:
        # No airgap found. This must be an inner or outer part
        logger.info("=== no airgap found ===")

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
                                  plt=p,  # plot
                                  kind=name,
                                  is_inner=inner,
                                  is_outer=outer,
                                  mindist=mindist,
                                  symtol=symtol,
                                  sympart=sympart,
                                  show_plots=show_plots,
                                  rows=3,  # rows
                                  cols=2,  # cols
                                  num=3)   # start num
        machine.delete_tiny_elements(mindist)

        if da > 0.0 or dy > 0.0:
            if inner:
                r_out = da / 2.0
                r_in = dy / 2.0
            elif outer:
                r_out = dy / 2.0
                r_in = da / 2.0
            else:
                r_out = 0.0
                r_in = 0.0
            if machine.cut_is_possible(r_in, r_out):
                logger.debug("make a cut")
                machine = machine.cut(r_in, r_out)

        if part:
            if part[0] == 'stator':
                machine.geom.set_stator()
                machine.set_inner_or_outer(part[1])
                machine.search_stator_subregions(single=True)
                machine = build_machine_stator(machine,
                                               inner,
                                               mindist,
                                               p,
                                               EESM=EESM,
                                               single=True)
                params = create_femag_parameters_stator(machine,
                                                        part[1])
            else:
                machine.geom.set_rotor()
                machine.set_inner_or_outer(part[1])
                machine.search_rotor_subregions(EESM, single=True)
                machine = build_machine_rotor(machine,
                                              inner,
                                              mindist,
                                              p,
                                              EESM=EESM,
                                              single=True)

                params = create_femag_parameters_rotor(machine,
                                                       part[1])
        else:
            machine.search_subregions(EESM, single=True)

        machine.create_inner_corner_areas()

        logger.info("***** END of work: %s *****", basename)

        if machine.geom.is_rotor():
            title = "Rotor"
        else:
            title = "Stator"

        if show_plots:
            p.render_elements(machine.geom, Shape,
                              draw_inside=True, title=title,
                              rows=3, cols=2, num=5, show=False,
                              with_corners=False,
                              with_nodes=False,
                              neighbors=False,
                              write_id=write_id,
                              draw_phi=True,
                              fill_areas=True)
        elif small_plots:
            #p.figure().suptitle(input_file.name, fontsize=16)
            p.render_elements(machine.geom, Shape,
                              draw_inside=True, title=title,
                              rows=1, cols=1, num=1, show=False,
                              with_corners=False,
                              with_nodes=False,
                              neighbors=False,
                              write_id=write_id,
                              draw_phi=True,
                              fill_areas=True)

        if show_plots or small_plots:
            if write_png:
                p.write_plot(basename)
            else:
                p.show_plot()

        if show_areas:
            p.render_elements(machine.geom, Shape,
                              title=name,
                              show=True,
                              draw_inside=True,
                              neighbors=True,
                              fill_areas=True)
            p.render_areas(machine.geom,
                           title=name,
                           with_nodes=True,
                           single_view=True)

        if write_fsl:
            if machine.is_full():
                logger.warning("it's not possible to create fsl-file")
                return None

            mtype = 'EESM' if EESM else 'PMSM'
            fslrenderer = FslRenderer(basename, mtype)
            conv['fsl'] = fslrenderer.render(machine, inner, outer, standalone=True)

    if params is not None:
        conv.update(params)

    if write_fsl:
        logger.debug("Write fsl")
        if conv and conv['fsl']:
            with io.open(basename + '.fsl', 'w',
                         encoding='utf-8') as f:
                f.write('\n'.join(conv['fsl']))
        else:
            logger.warning("No fsl data available")

        if conv:
            if conv.get('fsl_rotor_single', None):
                with io.open(basename + '_ROTOR.fsl', 'w',
                             encoding='utf-8') as f:
                    f.write('\n'.join(conv['fsl_rotor_single']))
                del conv['fsl_rotor_single']
            if conv.get('fsl_stator_single', None):
                with io.open(basename + '_STATOR.fsl', 'w',
                             encoding='utf-8') as f:
                    f.write('\n'.join(conv['fsl_stator_single']))
                del conv['fsl_stator_single']

    conv['name'] = basename
    t = timer.stop("-- all done in %0.4f seconds --", info=True)
    journal.put('time_total', t)
    journal.set('success', True)
    journal.write_journal()
    return conv


def _create_rotor_parameters(machine):
    rotor = {
        'min_radius': machine.geom.min_radius,
        'max_radius': machine.geom.max_radius,
        'mags': machine.geom.magnets_minmax_list(),
        'fd_wnds': machine.geom.fd_windings_minmax_list()
    }
    shaft_min, shaft_max = machine.geom.shaft_minmax()
    if shaft_max > 0.0:
        rotor['shaft_min'] = shaft_min
        rotor['shaft_max'] = shaft_max
    if shaft_max > rotor['min_radius']:
        rotor['min_radius'] = shaft_max
    return rotor


def _create_stator_parameters(machine):
    stator = {
        'min_radius': machine.geom.min_radius,
        'max_radius': machine.geom.max_radius,
        'wnds': machine.geom.windings_minmax_list()
    }
    shaft_min, shaft_max = machine.geom.shaft_minmax()
    if shaft_max > 0.0:
        stator['shaft_min'] = shaft_min
        stator['shaft_max'] = shaft_max
    if shaft_max > stator['min_radius']:
        stator['min_radius'] = shaft_max
    return stator


def create_femag_parameters(m_inner, m_outer, nodedist=1):
    if not (m_inner and m_outer):
        logger.warning("inner %s outer %s", m_inner, m_outer)
        return {}

    params = {}
    geom_inner = m_inner.geom
    geom_outer = m_outer.geom

    parts_inner = int(m_inner.get_symmetry_part())
    parts_outer = int(m_outer.get_symmetry_part())

    if m_inner.geom.is_rotor():
        geom_slots = geom_outer
        geom_poles = geom_inner
        num_slots = int(m_outer.get_symmetry_part())
        num_poles = int(m_inner.get_symmetry_part())
    else:
        geom_slots = geom_inner
        geom_poles = geom_outer
        num_slots = int(m_inner.get_symmetry_part())
        num_poles = int(m_outer.get_symmetry_part())

    slot_area = 0
    slot_area = geom_slots.area_size_of_type(TYPE_WINDINGS)
    num_sl_gen = int(geom_slots.get_symmetry_copies()+1)
    alfa_slot = geom_slots.get_alfa()
    alfa_pole = geom_poles.get_alfa()

    params['tot_num_slot'] = num_slots
    params['slot_area'] = slot_area
    params['num_sl_gen'] = num_sl_gen
    params['num_poles'] = num_poles
    params['nodedist'] = nodedist
    params['external_rotor'] = m_outer.geom.is_rotor()
    params['dy1'] = 2*geom_outer.max_radius
    params['da1'] = 2*geom_outer.min_radius
    params['da2'] = 2*geom_inner.max_radius
    params['dy2'] = 2*geom_inner.min_radius
    params['alfa_slot'] = alfa_slot
    params['alfa_pole'] = alfa_pole

    if m_inner.geom.is_rotor():
        params['rotor'] = _create_rotor_parameters(m_inner)
        params['stator'] = _create_stator_parameters(m_outer)
    else:
        params['rotor'] = _create_rotor_parameters(m_outer)
        params['stator'] = _create_stator_parameters(m_inner)

    if num_slots == 0 or num_poles == 0:
        if num_slots == 0:
            logger.warning("No slots found")
        if num_poles == 0:
            logger.warning("No poles found")
        logger.warning("Model not ready for femag")
        return {'error': "Model not ready for femag"}

    if not np.isclose(alfa_slot * num_slots,
                      alfa_pole * num_poles):
        logger.warning("slots and poles dont match")
        return {'error': "Model not ready for femag"}

    return params


def create_femag_parameters_stator(motor, position):
    params = {}
    num_slots = motor.get_num_slots()
    params['tot_num_slot'] = int(num_slots)
    if position == 'in':
        params['da2'] = 2*motor.geom.max_radius
        params['dy2'] = 2*motor.geom.min_radius
    else:
        params['dy1'] = 2*motor.geom.max_radius
        params['da1'] = 2*motor.geom.min_radius
    params['slot_area'] = motor.slot_area()
    params['stator'] = _create_stator_parameters(motor)
    params['machine'] = motor
    return params


def create_femag_parameters_rotor(motor, position):
    params = {}
    num_poles = motor.get_num_poles()
    params['num_poles'] = int(num_poles)
    if position == 'in':
        params['da2'] = 2*motor.geom.max_radius
        params['dy2'] = 2*motor.geom.min_radius
    else:
        params['dy1'] = 2*motor.geom.max_radius
        params['da1'] = 2*motor.geom.min_radius
    params['slot_area'] = motor.slot_area()
    params['rotor'] = _create_rotor_parameters(motor)
    params['machine'] = motor
    return params
