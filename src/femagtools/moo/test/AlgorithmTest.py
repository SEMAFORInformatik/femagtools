#!/usr/bin/env python
#
import unittest
import sys
sys.path.insert(0, '../..')
from moo.population import Population
from moo.problem import Problem
from moo.algorithm import Nsga2

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl
import matplotlib.cm as cm


class FesProblem(Problem):
    def __init__(self):
        super(FesProblem, self).__init__(4, 0, 2)
        
    def objfun(self, x):
        D = len(x)
        return (sum([abs(x[i]-np.exp(((i+1.)/D) ** 2)/3) ** 0.5
                    for i in range(D)]),
                sum([(x[i]-0.5*np.cos(10*(i+1)*np.pi/D)-0.5) ** 2
                    for i in range(D)]))


def surface():
    "show a nice 3d plot just for demo purposes"
    prob = FesProblem()
    steps = (15, 15)
    lower = (0.0, 0.0)
    upper = (1.0, 1.0)
    domain = [np.linspace(l, u, s)
              for s, l, u in zip(steps, lower, upper)]

    X = domain[0]
    Y = domain[1]

    L = [len(d) for d in domain]
    s = []
    e = 1
    for d in domain:
        L = L[1:]
        LS = np.prod(L)
        s.append(np.repeat(d.tolist()*LS, e))
        if len(L) > 0:
            e = e*L[0]

    z = np.array([prob.objfun(x) for x in np.array(s).T]).T[1]
    z = np.reshape(z, steps)

    fig = plt.figure()
    ax = mpl.Axes3D(fig)
    surf = ax.plot_surface(X, Y, z,
                           rstride=1,
                           cstride=1,
                           cmap=cm.jet,
                           linewidth=0)
    plt.show()


class AlgorithmTest(unittest.TestCase):
    def test_evolve(self):
        prob = FesProblem()
        pop = Population(prob, 20)
        pop.eval()

        F = np.array([ind.cur_f for ind in pop.individuals]).T
        plt.scatter(F[0], F[1], c="b", s=50, alpha=0.5)
        F = np.array([ind.cur_f
                      for ind in pop.individuals if ind.rank == 0]).T
        plt.scatter(F[0], F[1], c="g", s=100, alpha=0.5)
    
        algo = Nsga2()
        for i in range(25):
            newpop = algo.evolve(pop)
            newpop.eval()
            pop.merge(newpop)
            
        F = np.array([ind.cur_f for ind in pop.individuals]).T
        plt.scatter(F[0], F[1], c="r", s=200, alpha=0.5)
    
        plt.xlabel("$f^{(1)}$")
        plt.ylabel("$f^{(2)}$")
        plt.savefig('moo.png')

        self.assertEqual(pop.compute_norm_dist(), 1.0)
        self.assertEqual(pop.dom_count, [0]*pop.size())
        self.assertEqual(pop.pareto_rank, [0]*pop.size())

if __name__ == '__main__':
  unittest.main()
