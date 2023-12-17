# -*- coding: utf-8 -*-
"""
    femagtools.hxy
    ~~~~~~~~~~~~~~

    Reading HXY files (EXPERIMENTAL)
    Build cluster of magnets
"""

import numpy as np
import logging

# K-means clustering
# https://dev.to/sajal2692/coding-k-means-clustering-using-python-and-numpy-fg1
# Sajal Sharma
# with adaptions by Ronald Tanner
def initialize_random_centroids(K, X):
    """Initializes and returns k random centroids"""
    m, n = np.shape(X)
    # a centroid should be of shape (1, n), so the centroids array will be of shape (K, n)
    centroids = np.empty((K, n))
    for i in range(K):
        # pick a random data point from X as the centroid
        centroids[i] = X[np.random.choice(range(m))]
    return centroids


def initialize_centroids(K, X):
    """Initializes and returns k centroids"""
    c = np.mean(X, axis=0)
    if K < 2:
        return [c]
    alfa = np.arctan2(c[1], c[0])  # angle of center axis
    T = np.array([[np.cos(alfa), -np.sin(alfa)],
                  [np.sin(alfa), np.cos(alfa)]])
    p = (np.asarray(X)-c).dot(T) # rotate X
    d = np.linalg.norm(p, axis=1) # distance to center
    a = np.arctan2(p[:, 1], p[:, 0]) # angle to center
    h, b = np.histogram(a, bins=K)
    return [np.mean(X[(a>r[0]) & (a<r[1])], axis=0)
            for r in np.array((b[:-1], b[1:])).T]


def closest_centroid(x, centroids, K):
    """Finds and returns the index of the closest centroid for a given vector x"""
    distances = np.empty(K)
    for i in range(K):
        distances[i] = np.linalg.norm(centroids[i] - x)
    return np.argmin(distances) # return the index of the lowest distance


def create_clusters(centroids, K, X):
    """Returns an array of cluster indices for all the data samples"""
    m, _ = np.shape(X)
    cluster_idx = np.empty(m, dtype=int)
    for i in range(m):
        cluster_idx[i] = closest_centroid(X[i], centroids, K)
    return cluster_idx


def compute_means(cluster_idx, K, X):
    """Computes and returns the new centroids of the clusters"""
    _, n = np.shape(X)
    centroids = np.empty((K, n))
    for i in range(K):
        points = X[cluster_idx == i]  # gather points for the cluster i
        centroids[i] = np.mean(points, axis=0) # use axis=0 to compute means across points
    return centroids


def run_Kmeans(K, X, max_iterations=500):
    """Runs the K-means algorithm and computes the final clusters"""
    # initialize centroids
    centroids = initialize_centroids(K, X)
    # loop till max_iterations or convergance
    logging.debug(f"initial centroids: {centroids}")
    for _ in range(max_iterations):
        # create clusters by assigning the samples to the closet centroids
        clusters = create_clusters(centroids, K, X)
        previous_centroids = centroids
        # compute means of the clusters and assign to centroids
        centroids = compute_means(clusters, K, X)
        # if the new_centroids are the same as the old centroids, return clusters
        diff = previous_centroids - centroids
        if not diff.any():
            return clusters
    return clusters

"""references properties i to magnet k"""
class point():
    def __init__(self, index, k, coord):
        self.index = index
        self.coord = coord
        self.k = k

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
    """read hxy file and return values grouped to magnets
    returns:
      list of m=num_magnets dicts
        pos: list of n positions in degree
        e: n lists of center coordinates of elements in m
        hxy: n lists of field strengths in kA/m
        bxy: n lists of flux densities in T
        mxy: n lists of magnetization in T
        havg: average of field strengths kA/m
        hmax: maximum of field strengths kA/m
    """
    hxy = []
    with open(filename, encoding='latin1', errors='ignore') as f:
        n = 0
        k = 0
        for s in readSections(f):
            pos = float(s[0].split()[-1])
            num = np.array([[float(x) for x in l.split()] for l in s[5:] if l])
            hxy.append({'pos': pos, 'e': num[:, :2], 'hxy': num[:, 2:4],
                        'bxy': num[:, 4:6], 'mxy':num[:, 6:]})
            logging.info("HXY Section %d: pos %f shape %s", n, pos, num.shape)
            n += 1
            if k == 0:
                k = num.shape[1]
        K = num_magnets
        y_preds = run_Kmeans(num_magnets, np.array(hxy[0]['e']))
        logging.info("Kmeans: %s",
                     [y_preds[y_preds==k].shape[0] for k in range(K)])
        points = [point(i, k, hxy[0]['e'][i]) for i, k in enumerate(y_preds)]
        # move values to magnets:
        magnets = [{'e': [[h['e'][p.index]*1e-3
                           for p in points if p.k == k]
                    for h in hxy]}
                   for k in range(K)]

        for k, mag in enumerate(magnets):
            mag['pos'] = [h['pos'] for h in hxy]
            for mk in ['hxy', 'bxy', 'mxy']:
                mag[mk] = [[h[mk][p.index] for p in points if p.k == k]
                           for h in hxy]
            hxyabs = [np.linalg.norm(hxy, axis=1)
                      for hxy in mag['hxy']]
            mag['havg'] = [np.mean(a) for a in hxyabs]
            mag['hmax'] = [np.max(a) for a in hxyabs]

    return magnets

if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(message)s')
    magnets = read(sys.argv[1], int(sys.argv[2]))
    for m in magnets:
        print(f"{len(m['e'][0])}: Havg {m['havg']} Hmax {m['hmax']}")

    fig, ax = plt.subplots()
    for m in magnets:
        b = np.array(m['e'][0])
        ax.scatter(*b.T)
    ax.set_aspect('equal')
    ax.autoscale(enable=True)
    ax.axis('off')
    plt.show()
