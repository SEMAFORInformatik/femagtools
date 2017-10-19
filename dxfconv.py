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
import io

logger = logging.getLogger(__name__)


def usage(name):
    print("Usage: ", name,
          " [-h] [--help]", file=sys.stderr)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format='%(message)s')

    argparser = argparse.ArgumentParser(
        description='Process DXF file and create a plot or FSL file.')
    argparser.add_argument('dxfile', help='name of DXF file')
    argparser.add_argument('-f', help='create fsl', action="store_true")
    argparser.add_argument('-r', help='reshape based on symmetry detection',
                           action="store_true")
    args = argparser.parse_args()

    layers = ()
    
    incl_bnd = True
    
    pickdist = 1e-3
    basename = os.path.basename(args.dxfile).split('.')[0]
    logger.info("start reading %s", basename)

    basegeom = dg.Geometry(dg.dxfshapes(args.dxfile, layers=layers),
                           pickdist=pickdist)
    logger.info("total elements %s", len(basegeom.g.edges()))

    p = dr.PlotRenderer()

    motor_base = basegeom.get_motor()
    print("===== Original (nodes) =====")
    p.render_elements(basegeom, dg.Shape, with_nodes=True)
    if not motor_base.is_full():
        print("===== Original (hull) =====")
        p.render_elements(basegeom, dg.Shape, with_hull=True)
        print("===== Original (corners) =====")
        p.render_elements(basegeom, dg.Shape, with_corners=True)

        motor_base.repair_hull(basegeom)
        print("===== Original (REPAIRED HULL) =====")
        p.render_elements(basegeom, dg.Shape, with_corners=True)

    if not motor_base.is_a_motor():
        print("Not a Motor!!")
    else:
        geom_half = None
        geom_quarter = None
        
        if motor_base.is_full():
            motor_base.move_to_middle(basegeom)
            geom_half = motor_base.copy(basegeom, 0.0, np.pi)
        elif motor_base.is_half():
            motor_base.move_to_middle(basegeom)
            geom_half = motor_base.copy(basegeom, 0.0, 2*np.pi)
            
        if geom_half != None:
            motor_half = geom_half.get_motor()
            assert(motor_half.is_half())
            
            if not motor_half.is_in_middle():
                motor_half.move_to_middle(geom_half)
            
            motor_half.airgap(geom_half)
            print("===== Half =====")
            p.render_elements(geom_half, dg.Shape)
            
            geom_quarter = motor_half.copy(geom_half,0.0, np.pi/2)
        else:
            if motor_base.is_quarter():
                print("===== Quarter =====")
                motor_base.move_to_middle(basegeom)
                geom_quarter = motor_base.copy(basegeom, 0.0, 2.0*np.pi )
                
            else:
                # Es ist nicht klar, wie das Motorenteil aussieht
                print("Es ist nicht klar, wie das Motorenteil aussieht")
                motor_base.set_center(basegeom, 0.0, 0.0)
                motor_base.set_radius(9999999)
                geom_quarter = motor_base.copy(basegeom, 0.0, 2.0*np.pi )
                geom_quarter.clear_schnitt()
                print("===== Copy =====")
                p.render_elements(geom_quarter, dg.Shape)         
                print("===== Areas =====")            
                p.render_areas(geom_quarter, with_nodes=True)
            
        motor_quarter = geom_quarter.get_motor()
        if not motor_quarter.is_in_middle():
            motor_quarter.move_to_middle(geom_quarter)
        motor_quarter.rotate_to(geom_quarter, 0.0)
                        
        if motor_quarter.part_of_circle() == 0:
            print("Teil ist kein Teiler eines Kreises")
            p.render_elements(geom_quarter, dg.Shape)
            
        else:
            motor_quarter.airgap(geom_quarter)
            motor_quarter.repair_hull(geom_quarter)
            geom_quarter.clear_schnitt()
            print("===== Quarter =====")
            p.render_elements(geom_quarter, dg.Shape)
            
            if motor_quarter.has_airgap():
                geom_inner = motor_quarter.copy( geom_quarter, 0.0, 2*np.pi, True, True)
                motor_inner = geom_inner.get_motor()
                motor_inner.find_symmetrie(geom_inner)
                print("===== Inner =====")
                p.render_elements(geom_inner, dg.Shape)
                
                geom_outer = motor_quarter.copy( geom_quarter, 0.0, 2*np.pi, True, False)
                motor_outer = geom_quarter.get_motor()
                motor_outer.find_symmetrie(geom_outer)
                print("===== Outer =====")
                p.render_elements(geom_outer, dg.Shape)
            
#        p.render(geom_half, draw_corners=False, draw_center=False, incl_bnd=True, draw_hull=True)
        logger.info("done")
