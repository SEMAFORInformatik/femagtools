# -*- coding: utf-8 -*-
"""
    femagtools.parstudy
    ~~~~~~~~~~~~~~~~~~~

    Parameter study simulation



"""
import scipy as sc
import logging
import glob
import os
import time
import shutil
import functools
import pathlib
import numpy as np
import femagtools
import femagtools.model
import femagtools.fsl
import femagtools.condor
import femagtools.moproblem
import femagtools.getset
from .femag import set_magnet_properties

logger = logging.getLogger(__name__)


def get_report(decision_vars, objective_vars, objectives, x):
    """returns a combined list of objective and design values
    with header"""
    y = np.reshape(np.asarray(objectives),
                   (np.shape(objectives)[0], np.shape(x)[0])).T
    c = np.arange(x.shape[0]).reshape(x.shape[0], 1)
    return [[d['label'] for d in decision_vars] +
            [o['label']
             for o in objective_vars] + ['Directory'],
            [d['name'] for d in decision_vars] +
            [o['name']
             for o in objective_vars]] + [
        z[:-1].tolist() + [int(z[-1])]
        for z in np.hstack((x, y, c))]


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


class ParameterStudy(object):
    """abstract base class for parameter variation calculation

    Args:
    workdir: pathname of work directory

    CAUTION: choose the workdir exclusively for this usage.
             do not share it with anything else.
    """

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, condMat=[], result_func=None,
                 repname='grid', cmd=None,
                 templatedirs=[]):  # tasktype='Task'):
        #self.tasktype = tasktype
        self.templatedirs = templatedirs
        self.result_func = result_func
        self.femag = femagtools.Femag(workdir=workdir, cmd=cmd,
                                      magnetizingCurves=magnetizingCurves,
                                      magnets=magnets, condMat=condMat)
        # rudimentary: gives the ability to stop a running parameter variation. thomas.maier/OSWALD
        self.stop = False
        self.reportdir = ''
        self.repname = repname  # prefix for report filename ..-report.csv
        """
        the "owner" of the Grid have to take care to terminate all running xfemag64 or wfemagw64
        processes after setting stop to True
        For example:

        def killFemagThreads():
            if sys.platform.startswith('linux'):
                os.system(
                    "kill $(ps aux | grep '[x]femag64' | awk '{print $2}')")
            else:
                os.system("taskkill /f /im wfemagw64.exe")

        thomas.maier/OSWALD
        """

    def set_report_directory(self, dirname):
        """saves the result files (BCH/BATCH) of every calculation
        into this directory. Throws ValueError if directory is not empty or
        FileNotFoundError if it does not exist.
        Args:
          dirname: name of report directory
        """
        if os.listdir(dirname):
            raise ValueError("directory {} is not empty".format(dirname))
        self.reportdir = dirname

    def setup_model(self, builder, model, recsin):
        """builds model in current workdir and returns its filenames"""
        # get and write mag curves
        mc_files = self.femag.copy_magnetizing_curves(model, recsin=recsin)

        if model.is_complete():
            logger.info("setup model in %s", self.femag.workdir)
            filename = 'femag.fsl'
            if 'wdgdef' in model.windings:
                model.windings['wdgfile'] = self.femag.create_wdg_def(
                    model)

            with open(os.path.join(self.femag.workdir, filename), 'w') as f:
                f.write('\n'.join(builder.create_model(
                    model, self.femag.magnets, self.femag.condMat) +
                    ['save_model("close")']))

            self.femag.run(filename)

        model_files = [os.path.join(self.femag.workdir, m)
                       for m in mc_files] + [
            f for sublist in [
                glob.glob(os.path.join(
                    self.femag.workdir,
                    model.name+e))
                for e in (
                    '_*.poc', '*.WID',
                    '.nc', '.[IA]*7')]
            for f in sublist]

        return model_files

    def __call__(self, opt, machine, simulation,
                 engine, bchMapper=None,
                 extra_files=[], num_samples=0):
        """calculate objective vars for all decision vars
        Args:
          opt: variation parameter dict (decision_vars, objective_vars)
          machine: parameter dict of machine
          simulation: parameter dict of simulation
          engine: calculation runner (MultiProc, Condor ..)
          bchMapper: bch result transformation function
          extra_files: list of additional input file names to be copied
          num_samples: number of samples (ingored with Grid sampling)
        """

        self.stop = False  # make sure the calculation will start. thomas.maier/OSWALD

        decision_vars = opt['decision_vars']
        objective_vars = opt.get('objective_vars', {})

        model = femagtools.model.MachineModel(machine)
        builder = femagtools.fsl.Builder(self.templatedirs)

        dvarnames, domain, par_range = self._get_names_and_range(
            decision_vars, num_samples)

        # check if this model needs to be modified
        immutable_model = len([d for d in dvarnames
                               if hasattr(model,
                                          d.split('.')[0])]) == 0

        if isinstance(self.femag.workdir, pathlib.Path):
            workdir = self.femag.workdir
        else:
            workdir = pathlib.Path(self.femag.workdir)
        for d in workdir.glob('[0-9]*'):
            shutil.rmtree(d)

        extra_result_files = []
        if simulation.get('airgap_induc', False):
            extra_result_files.append('bag.dat')
        simulation['arm_length'] = model.lfe
        simulation['lfe'] = model.lfe
        simulation['move_action'] = model.move_action
        simulation['phi_start'] = 0.0
        simulation['range_phi'] = 720/model.get('poles')
        simulation.update(model.windings)
        if 'pocfilename' not in simulation:
            simulation['pocfilename'] = f"{model.name}_{model.poles}p.poc"
        fea = femagtools.model.FeaModel(simulation)

        prob = femagtools.moproblem.FemagMoProblem(decision_vars,
                                                   objective_vars)

        if immutable_model:
            modelfiles = self.setup_model(builder, model, recsin=fea.recsin)
            logger.info("Files %s", modelfiles+extra_files)

        job = engine.create_job(self.femag.workdir)
        # for progress logger
        job.num_cur_steps = fea.get_num_cur_steps()

        f = []
        p = 1
        calcid = 0
        logger.debug(par_range)

        if hasattr(fea, 'poc'):
            fea.poc.pole_pitch = 2*360/model.get('poles')
            fea.pocfilename = fea.poc.filename()
        elapsedTime = 0
        self.bchmapper_data = []  # clear bch data
        # split x value (par_range) array in handy chunks:
        popsize = 0
        for population in baskets(par_range, opt.get('population_size',
                                                     len(par_range))):
            if self.stop:  # try to return the results so far. thomas.maier/OSWALD
                logger.info(
                    'stopping grid execution... returning results so far...')
                try:
                    shape = [len(objective_vars)] + [len(d)
                                                     for d in reversed(domain)]
                    logger.debug("BEFORE: f shape %s --> %s",
                                 np.shape(np.array(f).T), shape)
                    completed = int(functools.reduce(
                        (lambda x, y: x * y), [len(z) for z in domain]))
                    logger.debug("need {} in total".format(completed))
                    remaining = completed - int(np.shape(np.array(f).T)[1])
                    values = int(np.shape(np.array(f).T)[0])
                    logger.debug("going to append {} None values".format(
                        remaining))
                    f += remaining * [values * [np.nan]]
                    shape = [len(objective_vars)] + [len(d)
                                                     for d in reversed(domain)]
                    logger.debug("AFTER: f shape %s --> %s",
                                 np.shape(np.array(f).T), shape)
                    objectives = np.reshape(np.array(f).T, shape)
                    r = dict(f=objectives.tolist(),
                             x=domain)
                    return r
                except:
                    return {}
                    pass
            popsize = max(len(population), popsize)
            logger.info('........ %d / %d results: %s',
                        p, int(np.ceil(len(par_range)/popsize)),
                        np.shape(f))
            job.cleanup()
            for k, x in enumerate(population):
                task = job.add_task(self.result_func)
                for fn in extra_files:
                    task.add_file(fn)
                if immutable_model:
                    prob.prepare(x, [fea, self.femag.magnets])
                    for m in modelfiles:
                        task.add_file(m)
                else:
                    prob.prepare(x, [model, fea, self.femag.magnets])
                    logger.info("prepare %s", x)
                    for mc in self.femag.copy_magnetizing_curves(
                            model,
                            dir=task.directory,
                            recsin=fea.recsin):
                        task.add_file(mc)

                set_magnet_properties(model, fea, self.femag.magnets)
                task.add_file(
                    'femag.fsl',
                    builder.create_model(model, self.femag.magnets) +
                    builder.create_analysis(fea) +
                    ['save_model("close")'])
                if hasattr(fea, 'poc'):
                    task.add_file(fea.pocfilename,
                                  fea.poc.content())
                if hasattr(fea, 'stateofproblem'):
                    task.set_stateofproblem(fea.stateofproblem)

            tstart = time.time()
            status = engine.submit(extra_result_files)
            logger.info('Started %s', status)
            status = engine.join()
            tend = time.time()
            elapsedTime += (tend-tstart)
            logger.info("Elapsed time %d s Status %s",
                        (tend-tstart), status)
            for t in job.tasks:
                if t.status == 'C':
                    r = t.get_results()
                    # save result file if requested:
                    if self.reportdir:
                        repdir = os.path.join(self.reportdir,
                                              str(calcid))
                        os.makedirs(repdir)
                        try:
                            shutil.copy(glob.glob(os.path.join(
                                t.directory, r.filename)+'.B*CH')[0], repdir)
                        except (FileNotFoundError, AttributeError):
                            # must be a failure, copy all files
                            for ff in glob.glob(
                                    os.path.join(t.directory, '*')):
                                shutil.copy(ff, repdir)
                        calcid += 1
                    if isinstance(r, dict) and 'error' in r:
                        logger.warn("job %d failed: %s", k, r['error'])
                        if objective_vars:
                            f.append([float('nan')]*len(objective_vars))
                        else:
                            f.append(dict())
                    else:
                        if bchMapper:
                            bchData = bchMapper(r)
                            self.addBchMapperData(bchData)
                            prob.setResult(bchData)
                        elif isinstance(r, dict):
                            prob.setResult(femagtools.getset.GetterSetter(r))
                        else:
                            prob.setResult(r)
                        f.append(prob.objfun([]))
                else:
                    if objective_vars:
                        f.append([float('nan')]*len(objective_vars))
                    else:
                        f.append(dict())
            p += 1

        logger.info('Total elapsed time %d s ...... DONE', elapsedTime)

        try:
            if objective_vars:
                # Note results f must be transposed but may have different shapes
                objectives = list(zip(*f))
            else:
                objectives = f
            if self.reportdir:
                self._write_report(decision_vars, objective_vars,
                                   objectives, par_range)
            return dict(f=objectives,
                        x=domain)
        except ValueError as v:
            logger.error(v)
            return dict(f=f, x=domain)

    def addBchMapperData(self, bchData):
        self.bchmapper_data.append(bchData)

    def getBchMapperData(self):
        return self.bchmapper_data

    def _write_report(self, decision_vars, objective_vars, objectives, par_range):
        with open(os.path.join(self.reportdir, f'{self.repname}-report.csv'), 'w') as f:
            for line in get_report(decision_vars, objective_vars,
                                   objectives, par_range):
                f.write(';'.join([str(v) for v in line]))
                f.write('\n')


class List(ParameterStudy):
    """List Parameter variation calculation"""

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, condMat=[],
                 result_func=None, cmd=None):  # tasktype='Task'):
        super(self.__class__, self).__init__(
            workdir,
            magnetizingCurves, magnets,
            condMat=condMat, result_func=result_func,
            repname='list', cmd=cmd)

    def _get_names_and_range(self, dvars, num_samples):
        if 'list' in dvars:
            par_range = dvars['list']
            domain = [r for r in zip(*par_range)]
            dnames = dvars['columns']
        else:
            domain = [d['values'] for d in dvars]
            par_range = [r for r in zip(*domain)]
            dnames = [d['name'] for d in dvars]
        return dnames, domain, par_range


class LatinHypercube(ParameterStudy):
    """Latin Hypercube sampling parameter variation calculation"""

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, condMat=[],
                 result_func=None, cmd=None):  # tasktype='Task'):
        super(self.__class__, self).__init__(
            workdir,
            magnetizingCurves, magnets,
            condMat=[], result_func=result_func, repname='lhs', cmd=cmd)

    def _get_names_and_range(self, dvars, num_samples):
        dvarnames = [d['name'] for d in dvars]
        l_bounds = [d['bounds'][0] for d in dvars]
        u_bounds = [d['bounds'][1] for d in dvars]

        N = num_samples
        sampler = sc.stats.qmc.LatinHypercube(d=len(l_bounds), centered=True)
        sample = sampler.random(n=N)
        par_range = sc.stats.qmc.scale(sample, l_bounds, u_bounds)
        domain = par_range.T.tolist()
        return dvarnames, domain, par_range


class Sobol(ParameterStudy):
    """Sobol sampling parameter variation calculation"""

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, condMat=[],
                 result_func=None, cmd=None):  # tasktype='Task'):
        super(self.__class__, self).__init__(workdir,
                                             magnetizingCurves, magnets,
                                             condMat=condMat, result_func=result_func,
                                             repname='sobol', cmd=cmd)

    def _get_names_and_range(self, dvars, num_samples):
        dvarnames = [d['name'] for d in dvars]
        l_bounds = [d['bounds'][0] for d in dvars]
        u_bounds = [d['bounds'][1] for d in dvars]

        N = num_samples
        sampler = sc.stats.qmc.Sobol(d=len(l_bounds), scramble=False)
        sample = sampler.random_base2(m=round(np.log(N)/np.log(2)))
        par_range = sc.stats.qmc.scale(sample, l_bounds, u_bounds)
        domain = par_range.T.tolist()
        return dvarnames, domain, par_range


class Grid(ParameterStudy):
    """Grid Parameter variation calculation"""

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, condMat=[],
                 result_func=None, cmd=None):  # tasktype='Task'):
        super(self.__class__, self).__init__(
            workdir,
            magnetizingCurves, magnets, condMat=condMat,
            result_func=result_func, cmd=cmd)

    def __create_parameter_range(self, domain):
        """returns the transposed array of the combined domain values"""
        L = [len(d) for d in domain]
        LS = np.prod(L)
        s = []
        e = 1
        for d in domain:
            LS = LS//len(d)
            s.append(np.repeat(d*LS, e))
            e = e*L[0]
            L = L[1:]
        return np.array(s).T

    def _get_names_and_range(self, dvars, num_samples):
        dvarnames = [d['name'] for d in dvars]
        domain = [list(np.linspace(d['bounds'][0],
                                   d['bounds'][1],
                                   d['steps']))
                  for d in dvars]

        par_range = self.__create_parameter_range(domain)
        return dvarnames, [r for r in zip(*par_range)], par_range
