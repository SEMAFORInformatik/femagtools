# -*- coding: utf-8 -*-
"""manage a geometry with
    lines, circles and arcs built from DXF

  NOTE: This code is in highly experimental state.
        Use at your own risk.

  Author: Ronald Tanner
    Date: 2017/07/06
"""
from __future__ import print_function
import numpy as np
import logging
from .functions import distance

logger = logging.getLogger('femagtools.corner')


#############################
#           Corner          #
#############################

class Corner(object):
    def __init__(self, center, p):
        self.__p = p
        self.__dist = distance(center, p)
        self.__keep = False
        self.is_new_point = False

    def point(self):
        return self.__p

    def is_equal(self, p, rtol=1e-04, atol=1e-04):
        return (np.isclose(self.__p[0], p[0], rtol=rtol, atol=atol) and
                np.isclose(self.__p[1], p[1], rtol=rtol, atol=atol))

    def is_same_corner(self, c):
        return self.is_equal(c.__p)

    def set_keep_node(self):
        self.__keep = True

    def keep_node(self):
        return self.__keep

    def __lt__(self, c):
        if self.__dist < c.__dist:
            return 1
        else:
            return 0

    def __str__(self):
        return "Corner: p={}".format(self.__p)
