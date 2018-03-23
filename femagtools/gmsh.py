# -*- coding: utf-8 -*-
"""
    femagtools.gmsh
    ~~~~~~~~~~~~~~~~

    Convert FEMAG modelfile to gmsh


"""
import femagtools.isa7
import logging
import os
import sys

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()

    isa = femagtools.isa7.read(filename)

    basename = os.path.splitext(os.path.basename(filename))[0]
    with open(basename + '.msh', 'w') as f:
        f.write('\n'.join(isa.msh()))
