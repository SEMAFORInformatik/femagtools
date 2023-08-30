# -*- coding: utf-8 -*-
"""
    dxfsl
    ~~~~~

    renders geometry


"""
import numpy as np
import logging

logger = logging.getLogger(__name__)


#############################
#       DumpRenderer        #
#############################

class DumpRenderer(object):
    def __init__(self, name):
        self.model = name

    def circle(self, center, radius):
        self.content.append(u'Crcle({}, {}, {})'.format(center[0],
                                                        center[1],
                                                        radius))

    def arc(self, startangle, endangle, center, radius):
        p1 = (center[0] + radius*np.cos(startangle),
              center[1] + radius*np.sin(startangle))
        p2 = (center[0] + radius*np.cos(endangle),
              center[1] + radius*np.sin(endangle))
        self.content.append(
            u"Arc({}, {}, {}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1],
                center[0], center[1]))

    def line(self, p1, p2):
        self.content.append(
            u"Line({}, {}, {}, {})".format(
                p1[0], p1[1], p2[0], p2[1]))

    def render(self, geom):
        '''create file with elements'''
        incl_bnd = True
        self.content = []
        handled = set()

        # collect all areas
        for i, area in enumerate(geom.areas(incl_bnd)):
            self.content.append('-- {}'.format(i))
            for e in area:
                if e not in handled:
                    e.render(self)
                    handled.add(e)
        self.content.append('--')
        geom.remove_areas(incl_bnd)
        # and all remaining edges and circles
        for circle in geom.circles():
            circle.render(self)

        for e1, e2, attr in geom.g.edges(data=True):
            attr['object'].render(self)

        return self.content
