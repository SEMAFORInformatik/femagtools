#!/usr/bin/env python3
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
#import io

logger = logging.getLogger(__name__)


def usage(name):
    print("Usage: ", name,
          " [-h] [--help]", file=sys.stderr)

def write_fsl(motor, basename):
    model = dr.NewFslRenderer(basename)
    filename = basename + '_' + motor.kind + '.fsl'
    model.render(motor.geom, filename)

def symmetry_search(motor, kind, sym_tolerance, show_plots):
    motor.clear_cut_lines()
    if show_plots:
        print("===== {} =====".format(kind))
        p.render_elements(motor.geom, dg.Shape, neighbors=True)
#        p.render_areas(motor.geom, with_nodes=False, single_view=True)
        
    if not motor.find_symmetry(sym_tolerance):
        print("symmetry_search: Keine Symmetrie gefunden")
        motor_mirror = motor.get_symmetry_mirror()
        motor_slice = motor
    else:
        p.render_elements(motor.geom, dg.Shape)
        motor_slice = motor.get_symmetry_slice()
        if motor_slice == None:
            motor.kind = kind
            return motor

        motor_mirror = motor_slice.get_symmetry_mirror()
        
    if motor_mirror == None:
        motor_ok = motor_slice
    else:
        if show_plots:
            print("===== Mirror of {} =====".format(kind))
            p.render_elements(motor_mirror.mirror_geom, dg.Shape)
            
        motor_ok = motor_mirror
                
    if show_plots:
        print("===== Final Result of {} =====".format(kind))
        p.render_elements(motor_ok.geom, dg.Shape, with_corners=True)
#        p.render_areas(motor_ok.geom, single_view=True)

    motor_ok.kind = kind
    return motor_ok

def write_fsl_file(motor):        
    if args.fsl:
        write_fsl(motor, basename, basename+"_"+kind+'.fsl')

   
#############################
#            Main           #
#############################
   
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='%(message)s')

    argparser = argparse.ArgumentParser(
        description='Process DXF file and create a plot or FSL file.')
    argparser.add_argument('dxfile',
                           help='name of DXF file')
    argparser.add_argument('-f', '--fsl',
                           help='create fsl',
                           action="store_true")
    argparser.add_argument('-a', '--airgap',
                           help='correct airgap',
                           dest='airgap',
                           type=float,
                           default=0.0)
    argparser.add_argument('-t', '--symtol',
                           help='absolut tolerance to find symmetrys',
                           dest='sym_tolerance',
                           type=float,
                           default=0.0)
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

    args = argparser.parse_args()
    if args.airgap > 0.0:
        print("Airgap is set to {}".format(args.airgap))
        
    layers = ()
    
    incl_bnd = True
    
    pickdist = 1e-3
    basename = os.path.basename(args.dxfile).split('.')[0]
    logger.info("start reading %s", basename)

    basegeom = dg.Geometry(dg.dxfshapes(args.dxfile, layers=layers),
                           pickdist=pickdist,
                           split=args.split)
    logger.info("total elements %s", len(basegeom.g.edges()))

    p = dr.PlotRenderer()
    if args.view:
        p.render_elements(basegeom, dg.Shape)
        sys.exit(0)

    motor_base = basegeom.get_motor()
    if args.show_plots:
        print("===== Original (nodes) =====")
        p.render_elements(basegeom, dg.Shape, neighbors=True)
        p.render_elements(basegeom, dg.Shape, with_hull=True)

    if not motor_base.is_a_motor():
        print("it's Not a Motor!!")
        sys.exit(1)
        
    if not motor_base.is_full():
        motor_base.repair_hull()
        if args.show_plots:
            print("===== Original (REPAIRED HULL) =====")
            p.render_elements(basegeom, dg.Shape, with_corners=True)

    if motor_base.is_full() or motor_base.is_half() or motor_base.is_quarter():
        # Wir erstellen eine Kopie des Originals für die weiteren Arbeiten
        motor = motor_base.full_copy()
    else:
        # Es ist nicht klar, wie das Motorenteil aussieht
        motor_base.set_center(0.0, 0.0)
        motor_base.set_radius(9999999)
        motor = motor_base.full_copy()
        
    if motor.part_of_circle() == 0:
        print("Teil ist kein Teilsegment eines Kreises")
        sys.exit(1)
        
    motor.clear_cut_lines()
    motor.move_to_middle()
    if args.show_plots:
        print("===== Areas =====")            
#        p.render_areas(motor.geom, with_nodes=True)
        p.render_elements(motor.geom, dg.Shape, with_corners=False)
        
    motor.airgap(args.airgap, args.sym_tolerance)
    
    if motor.has_airgap():
        motor_inner = motor.copy(0.0, 2*np.pi, True, True)
        motor_inner = symmetry_search(motor_inner, "Inner", args.sym_tolerance, args.show_plots)

        motor_outer = motor.copy(0.0, 2*np.pi, True, False)
        motor_outer = symmetry_search(motor_outer, "Outer", args.sym_tolerance, args.show_plots)
        motor_inner.sync_with_counterpart(motor_outer)
        
        if args.fsl:
            write_fsl(motor_inner, basename)
            write_fsl(motor_outer, basename)
        
    else:
        motor = symmetry_search(motor, "No_Airgap", args.sym_tolerance, args.show_plots)
        if args.fsl:
            write_fsl(motor, basename)
        
        
    logger.info("done")
