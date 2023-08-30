# -*- coding: utf-8 -*-
"""
    femagtools.dakotaout
    ~~~~~~~~~~~~~~~~~~~~

    parse dakota out file (helper module)

"""
import re
import logging

_numPattern = re.compile(r'\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)')
logger = logging.getLogger(__name__)


def _readSections(f):
    """return list of dakota out sections

    sections are split by empty lines

    Args:
      param f (file) dakota file to be read

    Returns:
      list of sections
    """

    section = []
    for line in f:
        if line.strip():
            section.append(line.rstrip())
        else:
            yield section
            section = []
    yield section


def read_moat_analysis(l):
    nums = _numPattern.findall(l)
    if len(nums) > 3:
        return [float(n) for n in nums[1:]]


def read_sample_moment(s):
    moments = []
    for l in s[2:]:
        nums = _numPattern.findall(l)
        if len(nums) > 3:
            name = l.split()[0]
            nums = _numPattern.findall(l[len(name):])
            moments.append(dict(
                name=name,
                mean=float(nums[0]),
                stddev=float(nums[1]),
                skewness=float(nums[2]),
                kurtosis=float(nums[3])))
    return 'moments', moments


def read_confidence(s):
    conf = []
    keys = ['LowerCI_Mean', 'UpperCI_Mean',
            'LowerCI_StdDev', 'UpperCI_StdDev']
    for l in s[2:]:
        nums = _numPattern.findall(l)
        if len(nums) > 3:
            name = l.split()[0]
            nums = _numPattern.findall(l[len(name):])

            d = {k: float(n) for k, n in zip(
                keys, nums)}
            d['name'] = name
            conf.append(d)
    return 'conf95', conf


def read_correlation(s):
    d = dict()
    k = ''
    for title in [('Simple Correlation', 'corr'),
                  ('Partial Correlation', 'corrpartial'),
                  ('Simple Rank Correlation', 'rank'),
                  ('Partial Rank Correlation', 'rankpartial')]:
        if s[0].startswith(title[0]):
            k = title[1]
            break
    if k:
        for l in s[2:]:
            nums = _numPattern.findall(l)
            if nums:
                d[l.split()[0]] = [float(n) for n in nums]
    return k, d


def read_pdf(s):
    pdf = []
    d = dict()
    name = ''
    for l in s[1:]:
        nums = _numPattern.findall(l)
        if len(nums) > 1:
            d['lower'].append(float(nums[0]))
            d['upper'].append(float(nums[1]))
            d['density'].append(float(nums[2]))
        elif l.startswith('PDF'):
            if d:
                pdf.append(d)
            d = dict(name=l.split()[-1][:-1],
                     lower=[], upper=[], density=[])
    return 'pdf', pdf


def read_cdf(s):
    cdf = []
    d = dict()
    name = ''
    for l in s[1:]:
        nums = _numPattern.findall(l)
        if len(nums) > 1:
            d['level'].append(float(nums[0]))
            d['p'].append(float(nums[1]))
        elif l.startswith('Cumulative'):
            if d:
                cdf.append(d)
            d = dict(name=l.split()[-1][:-1],
                     level=[], p=[])
    return 'cdf', cdf


def read_optsol(s):
    optsol = dict(parameters=[], objectives=[])
    for l in s[1:]:
        vn = l.split()
        if len(vn) == 2:
            nums = _numPattern.findall(f' {vn[0]}')
            if len(nums) == 1:
                optsol['parameters'].append((vn[1], float(nums[0])))
        elif len(vn) == 1:
            nums = _numPattern.findall(f' {vn[0]}')
            if len(nums) == 1:
                optsol['objectives'].append(float(nums[0]))

    return 'optsol', optsol


def read_dakota_out(filename):
    sectitles = [
        ('Sample moment statistics', read_sample_moment),
        ('95% confidence intervals', read_confidence),
        ('Probability Density Function', read_pdf),
        ('Level mappings', read_cdf),
        ('Simple Correlation', read_correlation),
        ('Partial Correlation', read_correlation),
        ('Simple Rank Correlation', read_correlation),
        ('Partial Rank Correlation', read_confidence),
        ('<<<<< Best parameters', read_optsol)]

    recs = dict()
    optsol = []
    with open(filename) as fp:
        stat = False
        for l in fp:
            if l.startswith('<<<<< Function'):
                break
        if fp and l:
            samples = int(_numPattern.findall(l)[0])
            recs['samples'] = samples
            for s in _readSections(fp):
                for t in sectitles:
                    if s:
                        if s[0].startswith(t[0]):
                            logger.info(s[0][:-1])
                            if t[1]:
                                r = t[1](s)
                                if r[0]:
                                    if r[0] == 'optsol':
                                        optsol.append(r[1])
                                    else:
                                        recs[r[0]] = r[1]
    if optsol:
        import numpy as np
        parnames = set([p[0] for o in optsol for p in o['parameters']])
        x = np.array([[p[1] for p in o['parameters']] for o in optsol]).T
        f = np.array([[f for f in o['objectives']] for o in optsol]).T
        recs['opt'] = dict(parnames=list(parnames), x=x.tolist(), f=f.tolist())
    return recs


if __name__ == '__main__':
    import sys
    import logging.config
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    for a in sys.argv[1:]:
        print(read_dakota_out(a))
