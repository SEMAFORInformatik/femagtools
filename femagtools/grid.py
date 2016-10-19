# -*- coding: utf-8 -*-
"""
    femagtools.grid
    ~~~~~~~~~~~~~~~

    Parameter range calculation

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import logging
import glob
import os
import numpy as np
import femagtools
import femagtools.model
import femagtools.fsl
import femagtools.condor
import femagtools.moproblem

logger = logging.getLogger(__name__)


def baskets(items, basketsize=10):
    """generates balanced baskets from iterable, contiguous items"""
    num_items = len(items)
    num_baskets = max(1, num_items//basketsize)
    if num_items % basketsize and basketsize < num_items:
        num_baskets += 1
    step = num_items//num_baskets
    for i in range(0, num_items, step):
        yield items[i:i+step]


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


class Grid(object):
    """Parameter variation calculation"""
    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None):
        self.femag = femagtools.Femag(workdir,
                                      magnetizingCurves=magnetizingCurves,
                                      magnets=magnets)

    def setup_model(self, builder, model):
        """builds model in current workdir and returns its filenames"""
        # get and write mag curves
        mc_files = self.femag.copy_magnetizing_curves(model)
        
        filename = 'femag.fsl'
        logger.info("setup model in %s", self.femag.workdir)
        with open(os.path.join(self.femag.workdir, filename), 'w') as f:
            f.write('\n'.join(builder.create_model(model,
                                                   self.femag.magnets) +
                              ['save_model(close)']))

        self.femag.run(filename, options=['-b'])
        model_files = [os.path.join(self.femag.workdir, m) for m in mc_files] + \
                      glob.glob(os.path.join(self.femag.workdir,
                                             model.name+'_*.poc')) + \
                        glob.glob(os.path.join(self.femag.workdir,
                                               model.name+'*7'))
        
        logger.info("model %s created", model.name)
        return model_files
    
    def create_parameter_range(self, domain):
        """returns the transposed array of the combined domain values"""
        L = [len(d) for d in domain]
        LS = np.prod(L)
        s = []
        e = 1
        for d in domain:
            LS = LS//len(d)
            s.append(np.repeat(d.tolist()*LS, e))
            e = e*L[0]
            L = L[1:]
        return np.array(s).T

    def __call__(self, opt, pmMachine, operatingConditions,
                 engine, bchMapper=None):
        """calculate objective vars for all decision vars"""
        decision_vars = opt['decision_vars']
        objective_vars = opt.get('objective_vars', {})

        steps = [d.get('steps', 10) for d in decision_vars]
        logger.info('STEPS %s', str(steps))

        model = femagtools.model.MachineModel(pmMachine)
        # check if this model needs to be modified
        immutable_model = len([d for d in decision_vars
                               if hasattr(model,
                                          d['name'].split('.')[0])]) == 0
        operatingConditions['lfe'] = model.lfe
        operatingConditions.update(model.windings)
        fea = femagtools.model.FeaModel(operatingConditions)

        prob = femagtools.moproblem.FemagMoProblem(decision_vars,
                                                   objective_vars)

        job = engine.create_job(self.femag.workdir)
        builder = femagtools.fsl.Builder()

        # build x value array
        domain = [np.linspace(l, u, s)
                  for s, l, u in zip(steps, prob.lower, prob.upper)]

        par_range = self.create_parameter_range(domain)
        f = []
        p = 1
        logger.debug(par_range)

        if immutable_model:
            modelfiles = self.setup_model(builder, model)
            logger.info("Files %s", modelfiles)

        self.bchmapper_data = []  # clear bch data
        # split x value (par_range) array in handy chunks:
        for population in baskets(par_range, opt['population_size']):
            logger.info('........ %d / %d', p, len(par_range)//len(population))
            if not isinstance(engine, femagtools.condor.Engine):
                job.cleanup()
            for k, x in enumerate(population):
                task = job.add_task()
                if immutable_model:
                    prob.prepare(x, fea)
                    for m in modelfiles:
                        task.add_file(m)
                    task.add_file('femag.fsl',
                                  builder.create_open(model) +
                                  builder.create_common(model) +
                                  builder.create_analysis(fea))
                else:
                    try:
                        prob.prepare(x, model)
                    except:
                        prob.prepare(x, [model, fea])
                    logger.info("prepare %s", x)
                    for mc in self.femag.copy_magnetizing_curves(
                            model,
                            task.directory):
                        task.add_file(mc)
                    task.add_file('femag.fsl',
                                  builder.create_model(model,
                                                       self.femag.magnets) +
                                  builder.create_analysis(fea))

            status = engine.submit()
            logger.info('Started %s', status)
            if bchMapper and isinstance(engine, femagtools.condor.Engine):
                 return {}  # BatchCalc Mode
            status = engine.join()

            for t in job.tasks:
                if t.status == 'C':
                    r = t.get_results()
                    if bchMapper:  # Mode => collectBchData
                        self.addBchMapperData(bchMapper(r))
                    if isinstance(r, dict) and 'error' in r:
                        logger.warn("job %d failed: %s", k, r['error'])
                        f.append([float('nan')]*len(objective_vars))
                    else:
                        prob.setResult(r)

                        f.append(prob.objfun([]))
            p += 1

        logger.info('...... DONE')
        logger.debug("Result %s", str(f))

        shape = [len(objective_vars)] + [len(d) for d in reversed(domain)]
        objectives = np.reshape(np.array(f).T, shape)
        return dict(f=objectives.tolist(),
                    x=[d.tolist() for d in domain])

    def addBchMapperData(self, bchData):
        delObjVecAttr = {'flux_fft': ['a', 'b', 'voltage_perc', 'flux_perc'],
                         'flux': ['displunit'],
                         'torque_fft': ['a', 'b']}
        delObjAttr = {'machineData': ['Q', 'p_sim', 'plfe', 'plfe1_0',
                                      'plfe2_0', 'plmag_0', 'qs_sim'],
                      'flux': ['displunit'],
                      'torque_fft': ['a', 'b'],
                      'lossPar': ['thetaw']}

        # remove object vector attributes
        for attr in delObjVecAttr:
            for i in range(len(bchData[attr])):
                for a in delObjVecAttr[attr]:
                    try:
                        del bchData[attr][i][a]
                    except:
                        pass

        # remove object attributes
        for attr in delObjAttr:
            for a in delObjAttr[attr]:
                try:
                    del bchData[attr][a]
                except:
                    pass

        self.bchmapper_data.append(bchData)

    def getBchMapperData(self):
        return self.bchmapper_data
