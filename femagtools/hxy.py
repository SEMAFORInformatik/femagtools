# -*- coding: utf-8 -*-
"""
    femagtools.hxy
    ~~~~~~~~~~~~~~

    Reading HXY files

"""

import numpy as np

def readSections(f):
    section = []
    movepos = False
    for line in f:
        if line.startswith(' MOVE POSITION'):
            movepos = True
            if section:
                # skip empty lines
                i = 0
                try:
                    while not section[i]:
                        i = i+1
                except IndexError:
                    i = i-1
                yield section[i:]
                section = []
        if movepos:
            section.append(line.strip())
    yield section

def read(filename):
    hxy = []
    with open(filename, encoding='latin1', errors='ignore') as f:
        for s in readSections(f):
            pos = float(s[0].split()[-1])
            num = np.array([[float(x) for x in l.split()] for l in s[5:] if l])
            h = np.linalg.norm(num[:, 2:4], axis=1)
            havg = np.mean(h)
            hmax = np.max(h)
            hxy.append({'pos': pos, 'e': num[:, :2], 'hxy': num[:, 2:4],
                        'bxy': num[:, 4:6], 'mxy':num[:, 6:],
                        'havg': havg, 'hmax': hmax})
    return hxy
