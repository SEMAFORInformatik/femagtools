# -*- coding: utf-8 -*-
"""
    femagtools.convert
    ~~~~~~~~~~~~~~~~~~

    Convert FEMAG files to various unstructured grid file formats
"""
import logging
import numpy as np
from femagtools import isa7
from collections import defaultdict

logger = logging.getLogger('femagtools.convert')


output_filetypes = [
    "msh", "geo", "vtu"]


def _from_isa(isa, filename, target_format,
              extrude=0, layers=0, recombine=False):

    if not isa.FC_RADIUS:
        logger.warn("airgap radius is not set in source file")

    airgap_center_elements = []
    for e in isa.elements:
        outside = [np.sqrt(v.x**2 + v.y**2) > isa.FC_RADIUS
                   for v in e.vertices]
        if any(outside) and not all(outside):
            airgap_center_elements.append(e)

    airgap_center_vertices = [v for e in airgap_center_elements
                              for v in e.vertices]

    airgap_inner_elements = []
    airgap_outer_elements = []
    for e in isa.elements:
        if e in airgap_center_elements:
            continue
        for v in e.vertices:
            if v not in airgap_center_vertices:
                continue
            if np.sqrt(v.x**2 + v.y**2) > isa.FC_RADIUS:
                airgap_outer_elements.append(e)
                break
            else:
                airgap_inner_elements.append(e)
                break

    airgap_outer_vertices = [v for e in airgap_outer_elements
                              for v in e.vertices]

    airgap_lines = []
    for e in airgap_center_elements:
        ev = e.vertices
        for i, v1 in enumerate(ev):
            v2 = ev[i-1]
            if v1 in airgap_outer_vertices and \
               v2 in airgap_outer_vertices:
                airgap_lines.append((v1, v2))

    nodechain_links = defaultdict(list)
    for nc in isa.nodechains:
        nodechain_links[nc.node1].extend(nc.nodes)
        nodechain_links[nc.node2].extend(nc.nodes)
        if nc.nodemid is not None:
            nodechain_links[nc.nodemid].extend(nc.nodes)

    physical_lines = ["v_potential_0",
                      "v_potential_const",
                      "periodic_+",
                      "periodic_-",
                      "infinite_boundary",
                      "no_condition",
                      "Airgap"]

    physical_surfaces = sorted(set([sr.name
                                    for sr in isa.subregions]
                                   + ["Winding_{}_{}".format(w.key, pol)
                                      for w in isa.windings
                                      for pol in ("+", "-")]
                                   + ["Air_Inner",
                                      "Air_Outer",
                                      "Airgap_Inner",
                                      "Airgap_Outer",
                                      "PM1", "PM2",
                                      "PM3", "PM4"]))

    def physical_line(n1, n2):
        if (n1, n2) in airgap_lines or (n2, n1) in airgap_lines:
            return 7  # airgap
        if n1.bndcnd == n2.bndcnd:
            return boundary_condition(n1)
        if boundary_condition(n1) == 1:
            return boundary_condition(n2)
        return boundary_condition(n1)

    def boundary_condition(node):
        if node.bndcnd == 0:
            return 6  # no condition
        if node.bndcnd == 1:
            return 1  # vpot 0
        if node.bndcnd == 2:
            return 2  # vpot const
        if node.bndcnd == 3 or node.bndcnd == 6:
            return 4  # periodic -
        if node.bndcnd == 4 or node.bndcnd == 5:
            return 3  # periodic +
        if node.bndcnd == 8 or node.bndcnd == 9:
            return 1  # vpot 0

    def physical_surface(e):

        def surface_id(name):
            return physical_surfaces.index(name) + len(physical_lines) + 1

        if any(e.mag):
            if e.mag[0] > 0:
                if e.mag[1] > 0:
                    return surface_id("PM1")
                return surface_id("PM2")
            if e.mag[1] > 0:
                return surface_id("PM3")
            return surface_id("PM4")

        if e in airgap_inner_elements or e in airgap_center_elements:
            return surface_id("Airgap_Inner")

        if e in airgap_outer_elements:
            return surface_id("Airgap_Outer")

        sr_key = isa.superelements[e.se_key].sr_key
        if sr_key == -1:
            v = e.vertices[0]
            if np.sqrt(v.x**2 + v.y**2) > isa.FC_RADIUS:
                return surface_id("Air_Outer")
            return surface_id("Air_Inner")

        sr = isa.subregions[sr_key]
        if sr.wb_key != -1:
            wb = isa.subregions[sr.wb_key]
            if sr.curdir > 0:
                return surface_id("Winding_{}_-".format(wb.key))
            return surface_id("Winding_{}_+".format(wb.key))

        return surface_id(sr.name)

    def line_on_boundary(n1, n2):
        if n1.on_boundary() and n2.on_boundary():
            if n1 in nodechain_links.keys():
                return n2 in nodechain_links[n1]
            else:
                return False
        else:
            return (n1, n2) in airgap_lines or (n2, n1) in airgap_lines

    points = [[n.x, n.y, 0] for n in isa.nodes]
    vpot = [n.vpot[0] for n in isa.nodes]

    lines = []
    line_ids = []
    triangles = []
    quads = []
    physical_ids = dict(triangle=[], quad=[])
    geometrical_ids = dict(triangle=[], quad=[])
    b = dict(triangle=[], quad=[])
    h = dict(triangle=[], quad=[])
    perm = dict(triangle=[], quad=[])
    iron_losses = dict(triangle=[], quad=[])
    mag_losses = dict(triangle=[], quad=[])
    wdg_losses = dict(triangle=[], quad=[])
    for e in isa.elements:
        
        ev = e.vertices
        
        for i, v in enumerate(ev):
            v1, v2 = v, ev[i-1]
            if line_on_boundary(v1, v2):
                lines.append([v1.key - 1, v2.key - 1])
                line_ids.append(physical_line(v1, v2))

        if len(ev) == 3:
            triangles.append([n.key - 1 for n in ev])
            cell_type = "triangle"

        elif len(ev) == 4:
            quads.append([n.key - 1 for n in ev])
            cell_type = "quad"
            
        physical_ids[cell_type].append(physical_surface(e))
        geometrical_ids[cell_type].append(e.se_key)
        b[cell_type].append(e.induction())
        h[cell_type].append(e.demagnetization())
        perm[cell_type].append(e.permeability())
        iron_losses[cell_type].append(e.iron_loss_density())
        mag_losses[cell_type].append(e.mag_loss_density())
        wdg_losses[cell_type].append(e.wdg_loss_density())

    if target_format == "msh":
        import meshio
        points = np.array(points)

        cells = {"line": np.array(lines),
                 "triangle": np.array(triangles),
                 "quad": np.array(quads)}

        point_data = {"potential": np.array(vpot)}

        cell_data = {
            "line": {
                "gmsh:geometrical": np.array(line_ids),
                "gmsh:physical": np.array(line_ids),
                "b": np.array([(0, 0, 0) for l in lines]),
                "h": np.array([0 for l in lines]),
                "Rel. Permeability": np.array([0 for l in lines]),
                "Iron Loss Dens.": np.array([0 for l in lines]),
                "Mag. Loss Dens.": np.array([0 for l in lines]),
                "Wdg. Loss Dens.": np.array([0 for l in lines])
            },
            "triangle": {
                "gmsh:geometrical": np.array(geometrical_ids["triangle"]),
                "gmsh:physical": np.array(physical_ids["triangle"]),
                "b": np.array([i + (0,) for i in b["triangle"]]),
                "h": np.array(h["triangle"]),
                "Rel. Permeability": np.array(perm["triangle"]),
                "Iron Loss Dens.": np.array(iron_losses["triangle"]),
                "Mag. Loss Dens.": np.array(mag_losses["triangle"]),
                "Wdg. Loss Dens.": np.array(wdg_losses["triangle"])
            },
            "quad": {
                "gmsh:geometrical": np.array(geometrical_ids['quad']),
                "gmsh:physical": np.array(physical_ids['quad']),
                "b": np.array([i + (0,) for i in b['quad']]),
                "h": np.array(h['quad']),
                "Rel. Permeability": np.array(perm['quad']),
                "Iron Loss Dens.": np.array(iron_losses['quad']),
                "Mag. Loss Dens.": np.array(mag_losses['quad']),
                "Wdg. Loss Dens.": np.array(wdg_losses['quad'])
            }
        }

        field_data = {}
        for l in physical_lines:
            field_data[l] = np.array([physical_lines.index(l) + 1, 1])
        for s in physical_surfaces:
            field_data[s] = np.array([physical_surfaces.index(s) + 1
                                      + len(physical_lines), 2])
        meshio.write_points_cells(filename,
                                  points,
                                  cells,
                                  point_data,
                                  cell_data,
                                  field_data,
                                  "gmsh2-ascii")

    if target_format == "geo":
        import meshio
        geo = []
        nc_nodes = set([n for nc in isa.nodechains for n in nc.nodes])

        for n in isa.nodes:
            if n in nc_nodes:
                geo.append("Point({}) = {{{}, {}, {}}};".format(
                    n.key, n.x, n.y, 0))

        for nc in isa.nodechains:
            geo.append("Line({}) = {{{}}};".format(
                    nc.key, ", ".join([str(n.key) for n in nc.nodes])))
        used = set()
        for nc in isa.nodechains:
            n1, n2 = nc.nodes[0], nc.nodes[1]
            if line_on_boundary(n1, n2):
                id_ = physical_line(n1, n2)
                name = physical_lines[id_ - 1]
                if extrude:
                    geo.append(
                        "extrusion[] = Extrude {{0, 0, {}}} {{ Line{{{}}}; {}{}}};".format(
                            extrude,
                            nc.key,
                            "Layers{{{}}}; ".format(layers) if layers else "",
                            "Recombine; " if recombine else ""))

                    geo.append("Physical Surface('{}', {}) {} extrusion[1];".format(
                        name,
                        id_,
                        "+=" if name in used else "="))
                else:
                    geo.append("Physical Line('{}', {}) {} {{{}}};".format(
                        name,
                        id_,
                        "+=" if name in used else "=",
                        nc.key))
                used.add(name)

        for se in isa.superelements:
            geo.append("Line Loop({}) = {{{}}};".format(
                se.key,
                ", ".join([str(nc.key) for nc in se.nodechains])))
            geo.append("Plane Surface({0}) = {{{0}}};".format(se.key))
        used = set()
        for se in isa.superelements:
            id_ = physical_surface(se.elements[0]) - len(physical_lines)
            name = physical_surfaces[id_ - 1]
            if extrude:
                geo.append(
                    "extrusion[] = Extrude {{0, 0, {}}} {{ Surface{{{}}}; {}{}}};".format(
                        extrude,
                        se.key,
                        "Layers{{{}}}; ".format(layers) if layers else "",
                        "Recombine; " if recombine else ""))

                geo.append("Physical Surface('base', {}) {} {{{}}};".format(
                    len(physical_lines) + 1,
                    "=" if se.key == 1 else "+=",
                    se.key))

                geo.append("Physical Surface('top', {}) {} extrusion[0];".format(
                    len(physical_lines) + 2,
                    "=" if se.key == 1 else "+="))

                geo.append("Physical Volume('{}', {}) {} extrusion[1];".format(
                    name,
                    id_,
                    "+=" if name in used else "="))
            else:
                geo.append("Physical Surface('{}', {}) {} {{{}}};".format(
                    name,
                    id_,
                    "+=" if name in used else "=",
                    se.key))
            used.add(name)
        
        with open(filename, "w") as f:
            f.write("\n".join(geo))

    if target_format == "vtu":
        import meshio
        assert len(points) == len(vpot)
        assert len(lines) == len(line_ids)
        assert len(triangles) == len(physical_ids['triangle']) == len(geometrical_ids['triangle']) == len(b['triangle']) == len(h['triangle']) == len(perm['triangle'])
        assert len(quads) == len(physical_ids['quad']) == len(geometrical_ids['quad']) == len(b['quad']) == len(h['quad']) == len(perm['quad'])

        points = np.array(points)

        cells = {"line": np.array(lines),
                 "triangle": np.array(triangles),
                 "quad": np.array(quads)}

        point_data = {"potential": np.array(vpot)}

        cell_data = {
            "line": {
                "GeometryIds": np.array(line_ids),
                "PhysicalIds": np.array(line_ids),
                "b": np.array([(0, 0) for l in lines]),
                "Demagnetization": np.array([0 for l in lines]),
                "Rel. Permeability": np.array([0 for l in lines]),
                "Iron Loss Dens.": np.array([0 for l in lines]),
                "Mag. Loss Dens.": np.array([0 for l in lines]),
                "Wdg. Loss Dens.": np.array([0 for l in lines])
            },
            "triangle": {
                "GeometryIds": np.array(geometrical_ids['triangle']),
                "PhysicalIds": np.array(physical_ids['triangle']),
                "b": np.array(b['triangle']),
                "Demagnetization": np.array(h['triangle']),
                "Rel. Permeability": np.array(perm['triangle']),
                "Iron Loss Dens.": np.array(iron_losses['triangle']),
                "Mag. Loss Dens.": np.array(mag_losses['triangle']),
                "Wdg. Loss Dens.": np.array(wdg_losses['triangle'])
            },
            "quad": {
                "GeometryIds": np.array(geometrical_ids['quad']),
                "PhysicalIds": np.array(physical_ids['quad']),
                "b": np.array(b['quad']),
                "Demagnetization": np.array(h['quad']),
                "Rel. Permeability": np.array(perm['quad']),
                "Iron Loss Dens.": np.array(iron_losses['quad']),
                "Mag. Loss Dens.": np.array(mag_losses['quad']),
                "Wdg. Loss Dens.": np.array(wdg_losses['quad'])
            }
        }

        field_data = {}
        for l in physical_lines:
            field_data[l] = np.array([physical_lines.index(l) + 1, 1])
        for s in physical_surfaces:
            field_data[s] = np.array([physical_surfaces.index(s) + 1
                                      + len(physical_lines), 2])
        meshio.write_points_cells(filename,
                                  points,
                                  cells,
                                  point_data=point_data,
                                  cell_data=cell_data,
                                  field_data=field_data,
                                  file_format="vtu-binary")


def to_msh(source, filename, infile_type=None):
    """
    Convert a femag model to msh format.

    Arguments:
        source: instance of femagtools.isa7.Isa7 or name of I7/ISA7 file
        filename: name of converted file
        infile_type: format of source file
    """
    if isinstance(source, isa7.Isa7):
        _from_isa(source, filename, "msh")

    elif type(source) == str:
        if infile_type:
            file_ext = infile_type.lower()
        else:
            file_ext = source.split(".")[-1].lower()

        if file_ext in ["isa7", "i7"]:
            isa = isa7.read(source)
            _from_isa(isa, filename, "msh")
        else:
            raise ValueError(
                "cannot convert files of format {} to .msh".format(file_ext))
    else:
        raise ValueError("cannot convert {} to .msh".format(source))

    
def to_geo(source, filename, extrude=0, layers=0,
           recombine=False, infile_type=None):
    """
    Convert a femag model to geo format.

    Arguments:
        source: instance of isa7.Isa7 or name of an I7/ISA7 file
        filename: name of converted file
        extrude: extrude surfaces using a translation along the z-axis
        layers: number of layers to create when extruding
        recombine: when extruding, recombine triangles into quadrangles and 
                 tetraedra into to prisms, hexahedra or pyramids
        infile_type: format of source file
    """
    if isinstance(source, isa7.Isa7):
        _from_isa(source, filename, "geo", extrude, layers, recombine)
        
    elif type(source) == str:
        if infile_type:
            file_ext = infile_type.lower()
        else:
            file_ext = source.split(".")[-1].lower()
        
        if file_ext in ["isa7", "i7"]:
            isa = isa7.read(source)
            _from_isa(isa, filename, "geo", extrude, layers, recombine)
        else:
            raise ValueError(
                "cannot convert files of format {} to .geo".format(file_ext))
    else:
        raise ValueError("cannot convert {} to .geo".format(source))


def to_vtu(source, filename, infile_type=None):
    """
    Convert a femag model to vtu format.

    Arguments:
        source: instance of isa7.Isa7 or name of an I7/ISA7 file
        filename: name of converted file
        infile_type: format of source file
    """
    if isinstance(source, isa7.Isa7):
        _from_isa(source, filename, "vtu")

    elif type(source) == str:
        if infile_type:
            file_ext = infile_type.lower()
        else:
            file_ext = source.split(".")[-1].lower()

        if file_ext in ["isa7", "i7"]:
            isa = isa7.read(source)
            _from_isa(isa, filename, "vtu")
        else:
            raise ValueError(
                "cannot convert files of format {} to .vtu".format(file_ext))
    else:
        raise ValueError("cannot convert {} to .vtu".format(source))


def main(argv=None):
    # Parse command line arguments.
    parser = _get_parser()
    args = parser.parse_args(argv)

    if not args.output_format:
        args.output_format = args.outfile.split('.')[-1]
        
    if args.output_format == 'msh':
        to_msh(args.infile, args.outfile)
    elif args.output_format == 'geo':
        to_geo(args.infile, args.outfile,
               extrude=args.extrude, layers=args.layers,
               recombine=args.recombine)
    elif args.output_format == 'vtu':
        to_vtu(args.infile, args.outfile)
    else:
        raise ValueError(
                "unsupported output format {}".format(args.output_format))
    return


def _get_parser():
    """Parse input options."""
    import argparse
    import sys
    from .__init__ import __version__
    
    parser = argparse.ArgumentParser(description=("Convert to mesh formats."))

    parser.add_argument("infile", type=str, help="mesh file to be read from")

    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=output_filetypes,
        help="output file format",
        default=None,
    )

    parser.add_argument(
        "--extrude",
        "-x",
        type=float,
        help="extrusion length",
        default=None,
    )
    
    parser.add_argument(
        "--layers",
        "-l",
        type=int,
        help="number of layers",
        default=None,
    )

    parser.add_argument(
        "--recombine",
        "-r",
        action='store_true',
        help="recombine")

    parser.add_argument("outfile", type=str, help="mesh file to be written to")

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version="%(prog)s {}, Python {}".format(__version__, sys.version),
        help="display version information",
    )

    return parser


if __name__ == "__main__":
    main()
