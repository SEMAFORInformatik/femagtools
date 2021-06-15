"""
  Copyright (c) 2015 Semafor Informatik & Energie AG

"""
import sys
import numpy as np
import copy
import operator as op


class Individual:
    def __init__(self, n_dim, f_dim):
        self.cur_f = [0]*f_dim
        self.cur_x = [0]*n_dim
        self.cur_v = []
        self.cur_c = []

        self.rank = 0
        self.crowd_d = 0
        self.idx = 0

    def __str__(self):
        return "f {} x {} rank {} crowd_d {}".format(self.cur_f, self.cur_x,
                                                     self.rank, self.crowd_d)


class Population:

    def __init__(self, probl, size, seed=None):
        self.individuals = []
        self.problem = probl
        self.champion = None
        self.dom_count = []
        self.dom_list = []
        np.random.seed(seed)
        for s in range(size):
            self.append([np.random.uniform(lb, ub)
                         for lb, ub in zip(self.problem.lower,
                                           self.problem.upper)])
#        kmax=int(np.sqrt(size))+1
#        f=[]
#        for i in range(kmax):
#            for k in range(kmax):
#                f.append((i*0.8/kmax +0.1, k*0.8/kmax +0.1))
#        for k,i in enumerate(self.individuals):
#            i.cur_f = f[k]

    def size(self):
        return len(self.individuals)

    def copy(self):
        return copy.copy(self)

    def append(self, child):
        self.individuals.append(Individual(self.problem.dimension,
                                           self.problem.f_dim))
        self.individuals[-1].cur_x = child
        self.individuals[-1].idx = len(self.individuals)-1
        self.init_velocity()

    def eval(self):
        for k, i in enumerate(self.individuals):
            i.cur_f = self.problem.objfun(i.cur_x)
        self.update()

    def best_idx(self):
        return [i.idx for i in sorted(self.individuals,
                                      key=op.attrgetter('rank',
                                                        'crowd_d'))]

    def merge(self, pop):
        "sort by rank and crowding distance (proximity)"
        self.individuals += pop.individuals
        self.update()
        for i in self.individuals:
            if abs(i.crowd_d) > sys.float_info.epsilon:
                i.crowd_d = 1/i.crowd_d
            else:
                i.crowd_d = sys.float_info.max

        best = sorted(self.individuals,
                      key=op.attrgetter('rank',
                                        'crowd_d'))
        self.individuals = best[:pop.size()]
        self.update()

    def init_velocity(self):
        for i, j in enumerate(self.individuals[-1].cur_x):
            w = (self.problem.upper[i] - self.problem.lower[i])/2
            self.individuals[-1].cur_v.append(np.random.uniform(-w, w))

    def populate(self, x, f, s):
        """populate with decision and objective values multiplied by sign
        x: vector of decision values
        f: matrix of objective values
        s: vector of signs (1, -1)"""
        for k, i in enumerate(self.individuals):
            i.cur_x = x[k]
            i.cur_f = [v*sign
                       for v, sign in zip(f[:, k], s)]
        self.update()

    def get_ranked_decisions(self):
        px = dict()
        for i in self.individuals:
            k = i.rank
            if k in px:
                px[k].append(i.cur_x)
            else:
                px[k] = [i.cur_x]
        return px

    def get_ranked_objectives(self, s):
        po = dict()
        for i in self.individuals:
            cur_f = [v*sign
                     for v, sign in zip(i.cur_f, s)]
            k = i.rank
            if k in po:
                po[k].append(cur_f)
            else:
                po[k] = [cur_f]
        return po

    def update(self):
        size = len(self.individuals)
        self.dom_count = []
        self.dom_list = [[] for s in range(size)]
        self.champion = None
        for s in range(size):
            self.dom_count.append(0)
            self.update_dom(s)
            self.update_champion(s)
        self.update_pareto_information()

    def update_dom(self, n):
        """Loop over the population (j) and construct
        dom_list[n] and dom_count """

        for i, x in enumerate(self.individuals):
            if i != n:
                # check if individual in position i
                # dominates the one in position n
                if self.problem.compare_fc(x.cur_f, x.cur_c,
                                           self.individuals[n].cur_f,
                                           self.individuals[n].cur_c):
                    self.dom_count[n] += 1
                    self.dom_list[i].append(n)

    def update_champion(self, idx):
        if self.champion is None or \
            self.problem.compare_fc(self.individuals[idx].cur_f,
                                    self.individuals[idx].cur_c,
                                    self.champion['f'], self.champion['c']):
            self.champion = dict(x=self.individuals[idx].cur_x,
                                 f=self.individuals[idx].cur_f,
                                 c=self.individuals[idx].cur_c)

    def update_crowding(self, F):
        # sort along the fitness dimension
        for i in range(self.problem.f_dim):
            I = sorted(F, key=lambda k: self.individuals[k].cur_f[i],
                       reverse=True)
            self.individuals[I[0]].crowd_d = sys.float_info.max
            self.individuals[I[-1]].crowd_d = sys.float_info.max
            df = (self.individuals[I[-1]].cur_f[i] -
                  self.individuals[I[0]].cur_f[i])
            for j in range(1, len(F)-1):
                if abs(df) > sys.float_info.epsilon:
                    self.individuals[I[j]].crowd_d += (
                        self.individuals[I[j+1]].cur_f[i] -
                        self.individuals[I[j-1]].cur_f[i])/df

    def update_pareto_information(self):
        size = len(self.individuals)
        self.pareto_rank = [0]*size
        F = [i for i, c in enumerate(self.dom_count) if c == 0]
        irank = 1
        dom_count_copy = list(self.dom_count)
        while True:
            self.update_crowding(F)
            S = []
            for i, f in enumerate(F):
                for j, k in enumerate(self.dom_list[f]):
                    dom_count_copy[k] -= 1
                    if dom_count_copy[k] == 0:
                        S.append(k)
                        self.pareto_rank[k] = irank
                        self.individuals[k].rank = irank
            if S:
                F = list(S)
            else:
                return

            irank += 1

    def compute_pareto_fronts(self):
        self.update_pareto_information()
        retval = [[] for s in range(max(self.pareto_rank)+1)]
        for i, j in enumerate(self.individuals):
            retval[self.pareto_rank[i]].append(i)
        return retval

    def compute_ideal(self):
        return [min(x)
                for x in zip(*[i.cur_f
                               for i in self.individuals if i.rank == 0])]

    def compute_nadir(self):
        return [max(x)
                for x in zip(*[i.cur_f
                               for i in self.individuals if i.rank == 0])]

    def compute_worst(self):
        return [max(x)
                for x in zip(*[i.cur_f
                               for i in self.individuals])]

    def compute_norm_dist(self):
        zi = np.array(self.compute_ideal())
        zw = np.array(self.compute_worst())
        znad = np.array(self.compute_nadir())
        return np.sqrt(((znad-zi)**2).sum() / ((zw-zi)**2).sum())

    def plot_pareto_fronts(
            self,
            rgb=(
                0,
                0,
                0),
            comp=[
                0,
                1],
            symbol='o',
            size=6,
            fronts=[]):
        """
        Plots the population pareto front in a 2-D graph

        USAGE: pop.plot_pareto_front(comp = [0,1], rgb=(0,1,0))

        * comp: components of the fitness function to plot in the 2-D window
        * rgb: specify the color of the 1st front (use strong colors here)
        * symbol: marker for the individual
        * size: size of the markersymbol
        * fronts: list of fronts to be plotted (use [0] to only show the first)
        """
        from numpy import linspace
        import matplotlib.pyplot as plt

        if len(comp) != 2:
            raise ValueError(
                'Invalid components of the objective function selected for plot')

        p_dim = self.problem.f_dimension

        if p_dim == 1:
            raise ValueError(
                'Pareto fronts of a 1-dimensional problem cannot be plotted')

        if not all([c in range(0, p_dim) for c in comp]):
            raise ValueError(
                'You need to select valid components of the objective function')

        p_list = self.compute_pareto_fronts()
        if (len(fronts) > 0):
            n = len(p_list)
            consistent = [d < n for d in fronts]
            if consistent.count(False) > 0:
                raise ValueError(
                    'Check your fronts list, there seem to be not enough fronts')
            p_list = [p_list[idx] for idx in fronts]

        cl = list(zip(linspace(0.9 if rgb[0] else 0.1, 0.9, len(p_list)),
                      linspace(0.9 if rgb[1] else 0.1, 0.9, len(p_list)),
                      linspace(0.9 if rgb[2] else 0.1, 0.9, len(p_list))))

        for id_f, f in enumerate(p_list):
            for ind in f:
                plt.plot([ind.cur_f[comp[0]]],
                         [ind.cur_f[comp[1]]],
                         symbol,
                         color=cl[id_f], markersize=size)
            x = [ind.cur_f[comp[0]] for ind in f]
            y = [ind.cur_f[comp[1]] for ind in f]
            tmp = [(a, b) for a, b in zip(x, y)]
            tmp = sorted(tmp, key=lambda k: k[0])
            plt.step([c[0] for c in tmp], [c[1]
                     for c in tmp], color=cl[id_f], where='post')
        return plt.gca()
