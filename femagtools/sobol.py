# -*- coding: utf-8 -*-
"""
    femagtools.sobol
    ~~~~~~~~~~~~~~~~

    Parameter range calculation using Sobol sampling



"""
import logging
import femagtools.sampling
import scipy as sc
import numpy as np

logger = logging.getLogger(__name__)


class Sobol(femagtools.sampling.Sampling):
    """Sobol sampling parameter variation calculation"""

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, result_func=None):  # tasktype='Task'):
        super(self.__class__, self).__init__(workdir,
                                             magnetizingCurves, magnets, result_func)

    def _get_names_and_range(self, dvars, num_samples):
        dvarnames = [d['name'] for d in dvars]
        l_bounds = [d['bounds'][0] for d in dvars]
        u_bounds = [d['bounds'][1] for d in dvars]

        N = num_samples
        sampler = sc.stats.qmc.Sobol(d=len(l_bounds), scramble=False)
        sample = sampler.random_base2(m=round(np.log(N)/np.log(2)))
        par_range = sc.stats.qmc.scale(sample, l_bounds, u_bounds)
        domain = par_range
        return dvarnames, domain, par_range
