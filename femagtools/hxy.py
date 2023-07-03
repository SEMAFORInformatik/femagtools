# -*- coding: utf-8 -*-
"""
    femagtools.hxy
    ~~~~~~~~~~~~~~

    Reading HXY files (EXPERIMENTAL)

"""

import numpy as np
from collections import defaultdict

# K-means clustering
class point():
    def __init__(self, index, k, coord):
        self.index = index
        self.coord = coord
        self.k = k

def make_k_mapping(points):
    region = defaultdict(list)
    for p in points:
        region[p.k] = region[p.k] + [p.coord]
    return region

def calc_k_means(region):
    return [np.mean(region[k], axis=0) for k in region]

def update_k(points, means):
    for p in points:
        dists = [np.linalg.norm(m - p.coord) for m in means]
        p.k = np.argmin(dists)

def fit(points, epochs=10):
    for e in range(epochs):
        region = make_k_mapping(points)
        means = calc_k_means(region)
        update_k(points, means)
    return means, points

def evaluate(points):
    region = make_k_mapping(points)
    means = calc_k_means(region)
    dists = [np.linalg.norm(means[p.k]-p.coord) for p in points]
    return np.mean(dists)

def llf_(y, X, pr):
    # return maximized log likelihood
    nobs = float(X.shape[0])
    nobs2 = nobs / 2.0
    nobs = float(nobs)
    resid = y - pr
    ssr = np.sum((resid)**2)
    llf = -nobs2*np.log(2*np.pi) - nobs2*np.log(ssr / nobs) - nobs2
    return llf


def aic(y, X, pr, p):
    # return aic metric
    llf = llf_(y, X, pr)
    return -2*llf+2*p


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


def read(filename, num_magnets):
    """read hxy file and return values grouped to magnets"""
    hxy = []
    with open(filename, encoding='latin1', errors='ignore') as f:
        for s in readSections(f):
            pos = float(s[0].split()[-1])
            num = np.array([[float(x) for x in l.split()] for l in s[5:] if l])
            hxy.append({'pos': pos, 'e': num[:, :2], 'hxy': num[:, 2:4],
                        'bxy': num[:, 4:6], 'mxy':num[:, 6:]})
        K = num_magnets
        points = [point(i, np.random.randint(0,K), xy)
                  for i, xy in enumerate(hxy[0]['e'])]
        new_means, new_points = fit(points)
        # move values to magnets:
        magnets = [{'e': [p.coord for p in new_points if p.k == k],
                    'pos': [], 'hxy': [], 'bxy': [], 'mxy': []}
                   for k in range(K)]
        hkeys = ['hxy', 'bxy', 'mxy']
        for i, h in enumerate(hxy):  # all positions
            for mag in magnets:
                mag['pos'].append(h['pos'])
                m = [{k: [] for k in hkeys}
                     for kk in range(K)]
            for p in new_points:  # all elements
                for k in hkeys:
                    m[p.k][k].append(h[k][p.k])
            for mk, magk in zip(m, magnets):
                for k in hkeys:
                    magk[k].append(mk[k])
        for mag in magnets:
            for k in ['e'] + hkeys:
                mag[k] = np.array(mag[k])
            mag['havg'] = []
            mag['hmax'] = []
            for hpos in mag['hxy']:
                h = np.abs(np.linalg.norm(hpos, axis=1))
                mag['havg'].append(np.mean(h))
                mag['hmax'].append(np.max(h))

        # Note dimension of hkeys is (positions x elements x 2)

    return magnets
