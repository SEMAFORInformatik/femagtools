"""read a dxf file and create a plot or fsl file

"""
import os
from pathlib import Path
from femagtools.dxfsl.geom import Geometry
from femagtools.dxfsl.shape import Shape
from femagtools.dxfsl.fslrenderer import FslRenderer, agndst
from femagtools.dxfsl.plotrenderer import PlotRenderer
from femagtools.dxfsl.concat import Concatenation
from femagtools.dxfsl.functions import Timer, middle_angle
from femagtools.dxfsl.journal import Journal, getJournal
import logging
import logging.config
import numpy as np
import sys

logger = logging.getLogger(__name__)
journal = None


def plot_geom(doit, plt, geom, title="Plot", areas=True):
    if not doit:
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
                    plt,  # plotter
                    kind,
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

    if is_inner:
        machine.set_inner()
    elif is_outer:
        machine.set_outer()

    machine.clear_cut_lines()
    if show_plots and debug_mode:
        plt.render_elements(machine.geom, Shape,
                            neighbors=True, title=kind)

    if sympart > 0:
        if not machine.is_full():
            logger.error("force symmetry failed")
            sys.exit(1)
        machine_ok = machine.get_forced_symmetry(sympart)
        machine_ok.set_kind(kind)
        logger.info("*** End of symmetry search for %s ***", kind)
        return machine_ok

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
            machine.kind = kind
            logger.info(" - no slice extracted ?!?")
            return machine

        plot_geom(False,  # for developer
                  plt, machine_slice.geom,
                  title="Slice of {}".format(kind))

        # no third after slice
        machine_mirror = machine_slice.get_symmetry_mirror(no_third=True)

        if machine_mirror is not None:
            machine_mirror.set_kind(kind)
            logger.info("*** End of symmetry search for %s (symmetry and mirror) ***", kind)
            return machine_mirror

        machine_slice.set_kind(kind)
        logger.info("*** End of symmetry search for %s (symmetry) ***", kind)
        return machine_slice

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
            machine_next_mirror = machine_mirror.get_symmetry_mirror()

        machine_ok = machine_mirror

    if not machine_ok.is_startangle_zero():
        logger.debug("Rotate geometry to 0.0")
        machine_ok.rotate_to(0.0)
        machine_ok.set_alfa_and_corners()

    machine_ok.set_kind(kind)

    logger.info("*** End of symmetry search for %s ***", kind)
    return machine_ok


def build_machine_rotor(machine, inner, mindist, plt, single=False):
    logger.debug("Begin of build_machine_rotor")
    if machine.has_windings():
        logger.debug("do nothing here with windings in rotor")
        logger.debug("End of build_machine_rotor")
        return machine

    if machine.is_mirrored():
        logger.debug("Rotor is mirrored")
        machine_temp = machine.undo_mirror()
        machine_temp.delete_tiny_elements(mindist)
        machine_temp.geom.set_rotor()
        machine_temp.search_rotor_subregions(single=single)
    else:
        machine_temp = machine

    midangle = middle_angle(machine_temp.startangle,
                            machine_temp.endangle)

    plot_geom(False,  # for developer
              plt, machine_temp.geom,
              title="Inner Rotor check magnets {}".format(midangle))

    rebuild = False
    if machine_temp.geom.magnets_in_the_middle(midangle):
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
        machine_temp.geom.area_list = []

    if machine_temp.create_auxiliary_lines():
        logger.debug("Auxiliary Lines created: rebuild subregions")
        rebuild = True
    if rebuild:
        machine_temp.rebuild_subregions(single=single)

    machine_temp.geom.recalculate_magnet_orientation()
    if inner:
        machine_temp.create_inner_corner_areas()

    if not machine_temp.is_mirrored():
        machine_temp.create_boundery_nodes()

    plot_geom(False,  # for developer
              plt, machine_temp.geom,
              title="Final Rotor")
    logger.debug("End of build_machine_rotor")
    return machine_temp


def build_machine_stator(machine, inner, mindist, plt, single=False):
    logger.debug("Begin of build_machine_stator")
    if not machine.geom.is_stator():
        logger.debug("Rotor with windings")

    if machine.has_mirrored_windings():
        logger.debug("undo mirrored windings")
        machine_temp = machine.undo_mirror()
        machine_temp.delete_tiny_elements(mindist)
        machine_temp.geom.set_stator()
        machine_temp.search_stator_subregions(single=single)
        machine_temp.create_mirror_lines_outside_windings()
    else:
        machine_temp = machine
    if machine_temp.create_auxiliary_lines():
        machine_temp.rebuild_subregions(single=single)

    if inner:
        machine_temp.create_inner_corner_areas()

    if not machine_temp.is_mirrored():
        machine_temp.create_boundery_nodes()

    plot_geom(False,  # for developer
              plt, machine_temp.geom,
              title="Final Stator")
    logger.debug("End of build_machine_stator")
    return machine_temp


def convert(dxfile,
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
            write_png=False,
            write_id=False,
            full_model=False,
            debug_mode=False):
    layers = ()
    conv = {}

    input_file = Path(dxfile)
    if not input_file.is_file():
        logger.error("File %s is not available", input_file)
        sys.exit(1)

    basename = input_file.stem
    if part:
        logger.info("***** start processing %s (%s) *****", basename, part)
    else:
        logger.info("***** start processing %s *****", basename)
    timer = Timer(start_it=True)

    journal = getJournal(name='converter', aktiv=debug_mode)
    journal.get_journal(input_file.name)
    journal.put_filename(str(input_file.resolve()))
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
    if small_plots:
        show_plots = False

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

        # inner part
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
        machine_inner.check_and_correct_geom("Inner")
        machine_inner.delete_tiny_elements(mindist)

        # outer part
        machine_outer = symmetry_search(machine_outer,
                                        p,  # plot
                                        outer_name,
                                        is_outer=True,
                                        symtol=symtol,
                                        show_plots=show_plots,
                                        rows=3,  # rows
                                        cols=2,  # columns
                                        num=4)   # start num
        machine_outer.check_and_correct_geom("Outer")
        machine_outer.delete_tiny_elements(mindist)

        machine_inner.sync_with_counterpart(machine_outer)

        machine_inner.search_subregions()
        machine_outer.search_subregions()

        # Inner mirrored rotor
        if machine_inner.geom.is_rotor():
            machine_inner = build_machine_rotor(machine_inner,
                                                True,  # is inner
                                                mindist,
                                                p)

        # Outer mirrored rotor
        if machine_outer.geom.is_rotor():
            machine_outer = build_machine_rotor(machine_outer,
                                                False,  # is outer
                                                mindist,
                                                p)

        if machine_inner.geom.is_stator() or machine_inner.has_windings():
            machine_inner = build_machine_stator(machine_inner,
                                                 True,
                                                 mindist,
                                                 p)

        if machine_outer.geom.is_stator() or machine_outer.has_windings():
            machine_outer = build_machine_stator(machine_outer,
                                                 False,
                                                 mindist,
                                                 p)
        machine_inner.sync_with_counterpart(machine_outer)

        machine_inner.search_critical_elements(mindist)
        machine_outer.search_critical_elements(mindist)

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

            fslrenderer = FslRenderer(basename)
            inner = fslrenderer.render(machine_inner, inner=True)
            outer = fslrenderer.render(machine_outer, outer=True)
            if full_model:
                params['num_sl_gen'] = params.get('tot_num_slot', 0)
            params['agndst'] = agndst(params.get('da1'),
                                      params.get('da2'),
                                      params.get('tot_num_slot'),
                                      params.get('num_poles'), 1)

            if params['external_rotor']:
                conv['fsl_magnet'] = outer
                conv['fsl_stator'] = inner
            else:
                conv['fsl_magnet'] = inner
                conv['fsl_stator'] = outer

            conv['fsl'] = fslrenderer.render_main(
                machine,
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
                                  p,  # plot
                                  name,
                                  is_inner=inner,
                                  is_outer=outer,
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
                logger.info("make a cut")
                machine = machine.cut(r_in, r_out)

        if part:
            if part[0] == 'stator':
                machine.geom.set_stator()
                machine.set_inner_or_outer(part[1])
                machine.search_stator_subregions(single=True)
                machine = build_machine_stator(machine, inner, mindist, p, single=True)

                params = create_femag_parameters_stator(machine,
                                                        part[1])
            else:
                machine.geom.set_rotor()
                machine.set_inner_or_outer(part[1])
                machine.search_rotor_subregions(single=True)
                machine = build_machine_rotor(machine, inner, mindist, p, single=True)

                params = create_femag_parameters_rotor(machine,
                                                       part[1])
        else:
            machine.search_subregions(single=True)

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

            fslrenderer = FslRenderer(basename)
            conv['fsl'] = fslrenderer.render(machine, inner, outer)

    if params is not None:
        conv.update(params)
    conv['name'] = basename
    t = timer.stop("-- all done in %0.4f seconds --", info=True)
    journal.put('time_total', t)
    journal.set('success', True)
    return conv


def create_femag_parameters(m_inner, m_outer, nodedist=1):
    if not (m_inner and m_outer):
        logger.warning("inner %s outer %s", m_inner, m_outer)
        return {}

    params = {}
    geom_inner = m_inner.geom
    geom_outer = m_outer.geom

    parts_inner = int(m_inner.get_symmetry_part())
    parts_outer = int(m_outer.get_symmetry_part())

    if parts_inner > parts_outer:
        num_slots = int(parts_inner)
        num_poles = int(parts_outer)
        num_sl_gen = int(geom_inner.get_symmetry_copies()+1)
        alfa_slot = geom_inner.get_alfa()
        alfa_pole = geom_outer.get_alfa()
    else:
        num_slots = int(parts_outer)
        num_poles = int(parts_inner)
        num_sl_gen = int(geom_outer.get_symmetry_copies()+1)
        alfa_slot = geom_outer.get_alfa()
        alfa_pole = geom_inner.get_alfa()

    params['tot_num_slot'] = num_slots
    params['num_sl_gen'] = num_sl_gen
    params['num_poles'] = num_poles
    params['nodedist'] = nodedist
    params['external_rotor'] = parts_inner > parts_outer
    params['dy1'] = 2*geom_outer.max_radius
    params['da1'] = 2*geom_outer.min_radius
    params['da2'] = 2*geom_inner.max_radius
    params['dy2'] = 2*geom_inner.min_radius
    params['alfa_slot'] = alfa_slot
    params['alfa_pole'] = alfa_pole

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
    return params
