# -*- coding: utf-8 -*-
"""
    femagtools.erg
    ~~~~~~~~~~~~~~

    Reading ERG files



"""
import re
import numpy as np
import io


head_pattern = re.compile('([A-Za-z_0-9]+)')
unit_pattern = re.compile('\[([A-Za-z/0-9 ]+)\]')
num_pattern = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)\s*')


def read(filename):
    """read ERG file
    returns dict with array values
    """
    head = []
    units = []
    m = []
    with io.open(filename,
                 errors='ignore') as f:
        for l in f:
            if head and units:
                n = num_pattern.findall(l)
                if n:
                    m.append([float(x) for x in n])
            elif head:
                u = unit_pattern.findall(l)
                if u:
                    units = u
            elif l.find('|') > 0:
                h = head_pattern.findall(l)
                if h:
                    head = h

    m = np.array(m).T
    ncols = len(set(m[1]))
    i1 = np.reshape(m[0], (-1, ncols)).T[0]
    nrows = len(i1)

    res = {k: (np.reshape(x, (nrows, ncols)).T[::-1]).tolist()
           for k, x in zip(head[2:], m[2:])}
    res['i1'] = i1.tolist()
    res['beta'] = m[1][:ncols][::-1].tolist()
    return res
