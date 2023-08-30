# -*- coding: utf-8 -*-
"""
    femagtools.gmsh
    ~~~~~~~~~~~~~~~~

    Handle gmsh models


"""
import numpy as np
import logging

class Gmsh(object):
    def __init__(self, mshfile):
        import meshio
        self.m = meshio.read(mshfile)
        self.r = np.linalg.norm(self.m.points, axis=1)
        self.phi = [np.arctan2(p[1], p[0]) for p in self.m.points]

    def get_section_angles(self):
        return np.min(self.phi), np.max(self.phi)

    def get_section_radius(self):
        return np.min(self.r), np.max(self.r)

    def get_subregions(self):
        return self.m.field_data.keys()

    def get_points(self, srname):
        """return x,y coordinates of all points in subregion srname"""
        srid = self.m.field_data[srname][0]
        trids = [i for i, k in enumerate(
            self.m.cell_data['gmsh:physical'][0]) if k == srid]
        quids = [i for i, k in enumerate(
            self.m.cell_data['gmsh:physical'][1]) if k == srid]
        p = np.unique(np.concatenate((self.m.cells[0].data[trids].flatten(), 
                                      self.m.cells[1].data[quids].flatten())))
        return self.m.points[p]
    
    def get_location(self, srname):
        """return x,y coordinates of first element (triangle, quad) in subregion srname"""
        srid = self.m.field_data[srname][0]
        eids = [i for i, k in enumerate(self.m.cell_data['gmsh:physical'][0]) if k == srid]
        if eids:
            return (np.sum(self.m.points[self.m.cells[0].data[eids[0]]], axis=0)/3)[:2]        
        eids = [i for i, k in enumerate(self.m.cell_data['gmsh:physical'][1]) if k == srid]
        if eids:
            return (np.sum(self.m.points[self.m.cells[1].data[eids[0]]], axis=0)/4)[:2]        
        raise ValueError("subregion '{}' not found".format(srname))

    def get_corners(self, srname):
        """return corner points in counterclockwise order"""
        p = self.get_points(srname)
        corner = []
        # 1: lower left
        for p1 in p:
            x = []
            for p2 in p:
                if p1[0] > p2[0]: 
                    x.append(p2)
                    break
            if not x:
                corner.append(p1[:2])
        # 2: lower right
        for p1 in p:
            x = []
            for p2 in p:
                if p1[1] > p2[1]: 
                    x.append(p2)
                    break
            if not x:
                corner.append(p1[:2])
        # 3: upper right
        for p1 in p:
            x = []
            for p2 in p:
                if p1[0] < p2[0]: 
                    x.append(p2)
                    break
            if not x:
                corner.append(p1[:2])
        # 4: upper left
        for p1 in p:
            x = []
            for p2 in p:
                if p1[1] < p2[1]: 
                    x.append(p2)
                    break
            if not x:
                corner.append(p1[:2])
        return corner

    def get_axis_angle(self, srname):
        """returns angle of axis in degrees"""
        corners = self.get_corners(srname)
        logging.debug("Corners %s of '%s'", corners, srname)
        l = corners[3], corners[2]
        alfa = np.arctan2(l[0][1] - l[1][1],
                          l[0][0] - l[1][0])/np.pi*180
                          
        logging.debug("Line l %s angle %s", l, alfa)
        return alfa
    
if __name__ == "__main__":
    import femagtools.isa7
    import logging
    import os
    import sys
    
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
