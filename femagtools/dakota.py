# -*- coding: utf-8 -*-
"""
    femagtools.dakota
    ~~~~~~~~~~~~~~~~~

    Running Dakota Analysis


"""
import sys
import logging
import json
import importlib
import pathlib
import subprocess
import femagtools
import mako
import mako.lookup
import numpy as np

logger = logging.getLogger(__name__)


class Dakota(object):
    """abstract base class for dakota calculations"""

    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, result_func=None,
                 template_name=''):
        # self.tasktype = tasktype
        self.result_func = result_func
        self.femag = femagtools.Femag(workdir,
                                      magnetizingCurves=magnetizingCurves,
                                      magnets=magnets)
        self.template_name = template_name
        if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
            # lookup up files in pyinstaller bundle
            logging.debug("Frozen!")
            dirs = [pathlib.Path(sys._MEIPASS) / 'dakotatemplates',
                    pathlib.Path.cwd()]
        else:
            dirs = [pathlib.Path(__file__).parent / 'templates/dakota',
                    pathlib.Path.cwd()]

        print([str(d) for d in dirs])
        self.lookup = mako.lookup.TemplateLookup(
            directories=[str(d) for d in dirs],
            disable_unicode=False, input_encoding='utf-8',
            output_encoding='utf-8',
            default_filters=['decode.utf8'])

    def __render(self, study, templ):
        try:
            template = self.lookup.get_template(templ+'.mako')
            logger.debug('use file {}'.format(templ))
            return template.render_unicode(study=study).split('\n')

        except mako.exceptions.TopLevelLookupException:
            logger.error('File {} not found'.format(templ))
            raise FileNotFoundError(templ)

    def _write_engine_conf(self, engine):
        # env f'ENGINE_{k.upper()}={engine[k]}' for k in engine]))
        (self.femag.workdir / 'engine.conf').write_text(
            '\n'.join([f'{k}={engine[k]}' for k in engine])+'\n')

    def __call__(self, study, machine, simulation,
                 engine, num_samples=0):
        """prepare input files and run dakota"""
        self._write_engine_conf(engine)
        (self.femag.workdir / 'model.py').write_text(
            '\n'.join([
                f'magnetMat={self.femag.magnets}',
                f'machine={machine}',
                f'simulation={simulation}']))

        model = femagtools.MachineModel(machine)
        self.femag.copy_magnetizing_curves(model)

        args = ['dakota', 'dakota.in']
        outname = self.femag.workdir / 'dakota.out'
        errname = self.femag.workdir / 'dakota.err'

        (self.femag.workdir / 'dakota.in').write_text(
            '\n'.join(self.__render(study, self.template_name)))

        with open(outname, 'w') as out, open(errname, 'w') as err:
            logger.info('invoking %s', ' '.join(args))
            proc = subprocess.Popen(
                args,
                stdout=out, stderr=err, cwd=self.femag.workdir)

        ret = proc.wait()
        if ret:
            raise RuntimeError(errname.read_text())
        varnames = ([d['name'] for d in study['decision_vars']] +
                    [o['name'] for o in study['objective_vars']])
        data = np.loadtxt(self.femag.workdir / 'femag.dat',
                          skiprows=1,
                          usecols=range(
                              2, 2 + len(varnames))).T
        xlen = len(study['decision_vars'])
        return dict(x=data[:xlen].tolist(),
                    f=data[xlen:].tolist())


class PsuadeMoat(Dakota):
    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, result_func=None):
        super(self.__class__, self).__init__(workdir,
                                             magnetizingCurves, magnets,
                                             result_func,
                                             template_name='psuade_moat')
