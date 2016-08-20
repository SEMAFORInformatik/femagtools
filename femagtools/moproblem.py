"""
    femagtools.moproblem
    ~~~~~~~~~~~~~~~~~~~~

    Creating and managing multi-objective problems

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
from .moo.problem import Problem
import logging

logger = logging.getLogger(__name__)


class FemagMoProblem(Problem):
    """
    A multi-objective optimization problem with Femag.

    """

    def __init__(self, decision_vars, objective_vars):
        Problem.__init__(self, len(decision_vars), 0, len(objective_vars))
        lbounds, ubounds = zip(*[d['bounds'] for d in decision_vars])
        self.set_bounds(lbounds, ubounds)
        self.decision_vars = decision_vars
        logger.info("Cnt decision_vars: %d", len(decision_vars))
        logger.info("bounds l: %s  u: %s", lbounds, ubounds)
        self.objective_vars = objective_vars

    # prepare model
    def prepare(self, x, model):
        # list of dicts (model, feaModel)
        if isinstance(model, list):
            for o in model:
                try:
                    self.prepare(x, o)
                except:
                    pass  # ignore silently
            return

        # simple dict
        for d, v in zip(self.decision_vars, x):
            logger.info("Prepare: %s = %s",  d['name'], v)
            model.set_value(d['name'].split('.'), v)

    def setResult(self, result):
        self.result = result
        if 'torque' in result.dqPar:
            self.result.dqPar['torque'] = result.dqPar['torque'][-1]
        self.result.torque = result.torque[-1]
        try:
            self.result.machine['plmag'] = result.machine['plmag'][-1]
            self.result.machine['plfe1'] = result.machine['plfe1'][-1]
            self.result.machine['plfe2'] = result.machine['plfe2'][-1]
            self.result.machine['plfe'] = result.machine['plfe'][-1]
        except:
            pass  # ignore any key or index errors
            
    def objfun(self, x):
        for o in self.objective_vars:
            logger.debug("%d=====> %s", len(self.objective_vars), str(o))
        return [f[0] * f[1] if f[1] else None
                for f in [(o.get('sign', 1),
                           self.result.get(o['name'].split('.')))
                          for o in self.objective_vars]]

    # Add some output to __repr__
    def __str__(self):
        return "\n\tMulti-Objective optimization problem"
