""" 
   Multiobjective Optimize Algorithm
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nondominated Sorting genetic algorithm II (NSGA-II)

 * NSGA-II is a nondominated-sorting based multiobjective evolutionary algorithm.
 * It genererates offspring with crossover and mutation and select the next
 * generation according to nondominated-sorting and crowding distance comparison.
 *
 * The algorithm can be applied to continuous box-bounded optimization. The version for mixed integer
 * and constrained optimization is also planned.
 * 
 * @see Deb, K. and Pratap, A. and Agarwal, S. and Meyarivan, T., 
       "A fast and elitist multiobjective genetic algorithm: NSGA-II"
 """
import random
import numpy as np
from .population import Population
import logging

logger = logging.getLogger("algorithm")


class Nsga2:
    def __init__(self):
        self.cr = 0.95
        self.eta_c = 10
        self.m = 0.01
        self.eta_m = 50
        pass

    def tournament_selection(self, i, j, pop):
        if pop.individuals[i].rank < pop.individuals[j].rank:
            return i
        if pop.individuals[i].rank > pop.individuals[j].rank:
            return j
        if pop.individuals[i].crowd_d < pop.individuals[j].crowd_d:
            return i
        if pop.individuals[i].crowd_d > pop.individuals[j].crowd_d:
            return j
        return random.sample((i, j), 2)[0]

    def crossover(self, pids, pop):
        "return 2 decision vectors"
        logger.debug("crossover {}".format(pids))

        parents = [pop.individuals[i] for i in pids]
        lb = pop.problem.lower
        ub = pop.problem.upper
        child1 = []
        child2 = []
        # a simulated binary crossover SBX
        if np.random.random_sample() <= self.cr:
            for x1, x2, yl, yu in zip(
                    parents[0].cur_x, parents[1].cur_x, lb, ub):
                if np.random.random_sample() <= 0.5 and \
                   abs(x1 - x2) > 1e-14:
                        if x1 < x2:
                            y1 = x1
                            y2 = x2
                        else:
                            y1 = x2
                            y2 = x1

                        rnd = np.random.random_sample()
                        beta = 1.0 + (2.0*(y1-yl)/(y2-y1))
                        alpha = 2.0 - beta**(-(self.eta_c+1.0))
                        if rnd <= (1.0/alpha):
                            betaq = (rnd*alpha)**(1.0/(self.eta_c+1.0))
                        else:
                            betaq = (1.0/(2.0 -
                                          rnd*alpha))**(1.0/(self.eta_c+1.0))

                        c1 = 0.5*((y1+y2)-betaq*(y2-y1))

                        beta = 1.0 + (2.0*(yu-y2)/(y2-y1))
                        alpha = 2.0 - beta**(-(self.eta_c+1.0))
                        if rnd <= (1.0/alpha):
                            betaq = (rnd*alpha)**(1.0/(self.eta_c+1.0))
                        else:
                            betaq = (1.0/(2.0 -
                                          rnd*alpha))**(1.0/(self.eta_c+1.0))

                        c2 = 0.5*((y1+y2)+betaq*(y2-y1))

                        if c1<yl : c1=yl
                        if c2<yl : c2=yl
                        if c1>yu : c1=yu
                        if c2>yu : c2=yu

                        if np.random.random_sample() <= 0.5:
                            child1.append(c1)
                            child2.append(c2)
                        else:
                            child1.append(c2)
                            child2.append(c1)
                else:
                    child1.append(x1)
                    child2.append(x2)
        else:
            child1 = parents[0].cur_x
            child2 = parents[1].cur_x
        return (child1, child2)

    def mutate(self, child, pop):
        mutant = []
        lb = pop.problem.lower
        ub = pop.problem.upper
        for i, x in enumerate(child):
            y = x
            yl = lb[i]
            yu = ub[i]
            delta1 = (y-yl)/(yu-yl)
            delta2 = (yu-y)/(yu-yl)
            rnd = np.random.random_sample()
            mut_pow = 1.0/(self.eta_m+1.0)
            if rnd<=0.5:
                xy = 1.0 - delta1
                val = 2.0*rnd+(1.0-2.0*rnd)*(xy**(self.eta_m+1.0))
                deltaq =  val**mut_pow - 1.0
            else:
                xy = 1.0 - delta2
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(xy**(self.eta_m+1.0))
                deltaq = 1.0 - val**mut_pow
            y = y + deltaq*(yu-yl)
            if y<yl: y = yl
            if y>yu: y = yu
            mutant.append(y)
        return mutant
    
    def evolve(self, pop):
        shuffles = (random.sample(range(pop.size()), pop.size()),
                    random.sample(range(pop.size()), pop.size()))
        newpop = Population(pop.problem, 0)
        for i in range(0, pop.size(), 4):
            for s in shuffles:
                parents = (self.tournament_selection(s[i], s[i+1], pop),
                           self.tournament_selection(s[i+2], s[i+3], pop))
                children = self.crossover(parents, pop)
                mutants = (self.mutate(children[0], pop),
                           self.mutate(children[1], pop))
                logger.debug("%s --> %s", children[0], mutants[0])
                logger.debug("%s --> %s", children[1], mutants[1])
                newpop.append(mutants[0])
                newpop.append(mutants[1])
        return newpop
