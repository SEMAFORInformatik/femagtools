# -*- coding: utf-8 -*-
"""
    femagtools.opt
    ~~~~~~~~~~~~~~

    Manage multi-objective optimization with FEMAG



"""
import femagtools
import femagtools.fsl
import femagtools.moproblem
from .moo.algorithm import Nsga2
from .moo.population import Population

import logging

logger = logging.getLogger(__name__)


def log_pop(pop, ngen):
    objectives = pop.problem.objective_vars
    decisions = pop.problem.decision_vars
    log = [''.join(['Generation: {}\n'.format(ngen),
                    'rank '] +
                   ["{:10s}".format(o['name'].split('.')[-1])
                    for o in objectives] +
                   [" "] +
                   ["{:12s}".format(d['name'].split('.')[-1])
                    for d in decisions])]

    for i in pop.individuals:
        log.append(''.join(["{},".format(i.rank)] +
                           ["{:10.2f},".format(f) for f in i.cur_f] +
                           ["   "] +
                           ["{:10.4f},".format(x) for x in i.cur_x]))

    log.append(
        '--------------------------------------------------------------------')
    log.append("  znad: {}\n  zi: {}\n  zw: {}\n  Norm Dist: {}".format(
        pop.compute_nadir(),
        pop.compute_ideal(),
        pop.compute_worst(),
        pop.compute_norm_dist()))
    return log


class Optimizer(object):
    def __init__(self, workdir, magnetizingCurves, magnetMat):
        self.femag = femagtools.Femag(workdir,
                                      magnetizingCurves=magnetizingCurves,
                                      magnets=magnetMat)
        
    def _update_population(self, pop, engine):
        self.job.cleanup()

        for k, i in enumerate(pop.individuals):
            task = self.job.add_task()
            pop.problem.prepare(i.cur_x, self.model)
            for mc in self.femag.copy_magnetizing_curves(self.model,
                                                         task.directory):
                task.add_file(mc)
            task.add_file('femag.fsl',
                          self.builder.create(self.model, self.fea,
                                              self.femag.magnets))

        ntasks = engine.submit()
        status = engine.join()

        for t, i in zip(self.job.tasks, pop.individuals):
            if t.status == 'C':
                r = t.get_results()
                if isinstance(r, dict) and 'error' in r:
                    logger.warn("Task %s failed: %s", t.id, r['error'])
                else:
                    pop.problem.setResult(r)

                    i.cur_f = pop.problem.objfun([])
            else:
                logger.warn("Task %s failed with status %s", t.id, t.status)

        pop.update()

    def __call__(self, num_generations, opt, pmMachine,
                 operatingConditions, engine):
        return self.optimize(num_generations, opt, pmMachine,
                             operatingConditions, engine)
        
    def optimize(self, num_generations, opt, pmMachine,
                 operatingConditions, engine):
        """execute optimization"""
        decision_vars = opt['decision_vars']
        objective_vars = opt['objective_vars']
        population_size = opt['population_size']
                    
        problem = femagtools.moproblem.FemagMoProblem(decision_vars,
                                                      objective_vars)
        self.builder = femagtools.fsl.Builder()
        self.model = femagtools.model.MachineModel(pmMachine)
        self.fea = operatingConditions
        self.fea.update(self.model.windings)
        self.fea['lfe'] = self.model.lfe
        self.fea['move_action'] = self.model.move_action
        self.fea['pocfilename'] = (self.model.get('name') +
                                   '_' + str(self.model.get('poles')) +
                                   'p.poc')
        self.pop = Population(problem, population_size)
        
        algo = Nsga2()

        self.job = engine.create_job(self.femag.workdir)
        
        logger.info("Optimize x:%d f:%d generations:%d population size:%d",
                    len(self.pop.problem.decision_vars),
                    len(self.pop.problem.objective_vars),
                    num_generations,
                    self.pop.size())

        results = dict(rank=[], f=[], x=[])

        for i in range(num_generations):
            logger.info("Generation %d", i)
            if i > 0:
                newpop = algo.evolve(self.pop)
                self._update_population(newpop, engine)
                self.pop.merge(newpop)
            else:
                self._update_population(self.pop, engine)
            logger.info('\n'.join(log_pop(self.pop, i)))
        logger.info("finished")
        ft = []
        xt = []
        for i in self.pop.individuals:
            results['rank'].append(i.rank)
            ft.append(i.cur_f)
            xt.append(i.cur_x)
        objective_vars = self.pop.problem.objective_vars
        decision_vars = self.pop.problem.decision_vars
        results['f'] = [[s.get('sign', 1)*y for y in f]
                        for s, f in zip(objective_vars, zip(*ft))]
        results['x'] = list(zip(*xt))
        results['znad'] = [s.get('sign', 1)*v
                           for s, v in zip(objective_vars,
                                           self.pop.compute_nadir())]
        results['zi'] = [s.get('sign', 1)*v
                         for s, v in zip(objective_vars,
                                         self.pop.compute_ideal())]
        results['zw'] = [s.get('sign', 1)*v
                         for s, v in zip(objective_vars,
                                         self.pop.compute_worst())]
        results['dist'] = self.pop.compute_norm_dist()
        results['objective'] = [o['desc'] for o in objective_vars]
        results['decision'] = [d['desc'] for d in decision_vars]
        return results
    
    #print("\nChampion: {}\n".format(pop.champion['f']))
        #if flast != None:
        #    print("Fitness Comparison:")
        #    for f1, f2 in zip(pop.champion['f'], flast):
        #        print( "{:10.2f} {:10.2f}      {:10.2f}".format(f1, f2, f1-f2))
        #flast = list(pop.champion['f'])
        #print("")
#    except:
#        print "Unexpected error:", sys.exc_info()

#    print("L2 distance to the best decision vector:")
#    for best_decision in prob.best_x:
#        l2_norm = 0
#        for n in range(0, len(best_decision)):
#            l2_norm +=  (best_decision[n] - isl.population.champion.x[n]) ** 2
#        l2_norm = sqrt(l2_norm)
#        print(l2_norm)

#pf=pop.plot_pareto_fronts()
#savefig('pf0.png')
