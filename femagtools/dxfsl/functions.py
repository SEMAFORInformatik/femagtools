# -*- coding: utf-8 -*-
"""manage a geometry with
    lines, circles and arcs built from DXF

  NOTE: This code is in highly experimental state.
        Use at your own risk.

  Author: Ronald Tanner
    Date: 2017/07/06
"""
from __future__ import print_function
import logging
import numpy as np

logger = logging.getLogger('femagtools.functions')


def distance(p1, p2):
    assert(len(p1) > 1)
    assert(len(p2) > 1)
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
