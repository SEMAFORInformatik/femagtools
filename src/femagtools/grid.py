# -*- coding: utf-8 -*-
"""
    femagtools.grid
    ~~~~~~~~~~~~~~~

    Parameter range calculation using Grid sampling



"""
import logging
import femagtools.parstudy
import numpy as np

logger = logging.getLogger(__name__)


def create_parameter_range(domain):
    """returns the transposed array of the combined domain values"""
    L = [len(d) for d in domain]
    LS = np.prod(L)
    s = []
    e = 1
    for d in domain:
        LS = LS//len(d)
        s.append(np.repeat(d*LS, e))
        e = e*L[0]
        L = L[1:]
    return np.array(s).T


class Grid(femagtools.parstudy.ParameterStudy):
    """Grid Parameter variation calculation"""

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, condMat=[],
                 result_func=None):  # tasktype='Task'):
        super().__init__(workdir,
                         magnetizingCurves, magnets, condMat,
                         result_func)

    def _get_names_and_range(self, dvars, num_samples):
        if isinstance(dvars, dict):
            dvarnames = dvars['columns']
            par_range = dvars['list']
            domain = [r for r in zip(*par_range)]
        else:
            dvarnames = [d['name'] for d in dvars]
            domain = [list(np.linspace(d['bounds'][0],
                                       d['bounds'][1],
                                       d['steps'])) if 'bounds' in d else d['values']
                      for d in dvars]

            par_range = create_parameter_range(domain)
        return dvarnames, domain, par_range
