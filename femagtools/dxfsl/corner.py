# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.corner
  ~~~~~~~~~~~~~~~~~~~~~~~

  Authors: Ronald Tanner, beat Holm
"""
from __future__ import print_function
import numpy as np
import logging
from .functions import distance

logger = logging.getLogger('femagtools.corner')


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

    def __eq__(self, c):
        return (np.isclose(self.__p[0], c.__p[0], rtol=1e-04, atol=1e-04) and
                np.isclose(self.__p[1], c.__p[1], rtol=1e-04, atol=1e-04))

    def __lt__(self, c):
        if self.__dist < c.__dist:
            return 1
        else:
            return 0

    def __str__(self):
        return "Corner: p={}, keep={}".format(self.__p, self.__keep)
