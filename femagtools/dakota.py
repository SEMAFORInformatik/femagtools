# -*- coding: utf-8 -*-
"""
    femagtools.dakota
    ~~~~~~~~~~~~~~~~~

    Running Dakota Analysis


"""
import sys
import os
import logging
import copy
import pathlib
import subprocess
import femagtools
import mako
import mako.lookup
import numpy as np
import json
import asyncio
from asyncio.subprocess import PIPE
from .dakotaout import read_dakota_out

logger = logging.getLogger(__name__)

use_asyncio = True


@asyncio.coroutine
def read_stream_and_display(stream, display):
    """Read from stream line by line until EOF, display, and capture the lines.

    """
    output = []
    while True:
        line = yield from stream.readline()
        if not line:
            break
        output.append(line)
        display(line)  # assume it doesn't block
    return b''.join(output)


@asyncio.coroutine
def read_and_display(*cmd, cwd):
    """Capture cmd's stdout, stderr while displaying them as they arrive
    (line by line).

    """
    # start process
    process = yield from asyncio.create_subprocess_exec(*cmd, cwd=cwd,
                                                        stdout=PIPE, stderr=PIPE)

    def nodisplay(l):
        return
    # read child's stdout/stderr concurrently (capture and display)
    try:
        stdout, stderr = yield from asyncio.gather(
            read_stream_and_display(process.stdout, nodisplay),
            read_stream_and_display(process.stderr, sys.stderr.buffer.write))
    except Exception:
        process.kill()
        raise
    finally:
        # wait for the process to exit
        rc = yield from process.wait()
    return rc, stdout, stderr


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
                 engine, **kwargs):
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

        dstudy = copy.deepcopy(study)
        dstudy.update(kwargs)
        for o in dstudy['objective_vars']:
            if o.get('sign', 1) < 0:
                o['name'] = f'-{o["name"]}'

        (self.femag.workdir / 'dakota.in').write_text(
            '\n'.join(self.__render(dstudy, self.template_name)))

        logger.info('invoking %s with %s',
                    ' '.join(args), engine['module'])
        if use_asyncio:
            # run the event loop
            if os.name == 'nt':
                loop = asyncio.ProactorEventLoop()  # for subprocess' pipes on Windows
                asyncio.set_event_loop(loop)
            else:
                loop = asyncio.get_event_loop()
            ret, *output = loop.run_until_complete(
                read_and_display(*args, cwd=self.femag.workdir))
            loop.close()
            outname.write_text(output[0].decode())
            errname.write_text(output[1].decode())
        else:
            with open(outname, 'w') as out, open(errname, 'w') as err:
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
        result = dict(x=data[:xlen].tolist(),
                      f=data[xlen:].tolist())
        result.update(read_dakota_out(
            self.femag.workdir / 'dakota.out'))

        return result


class Sampling(Dakota):
    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, result_func=None):
        super(self.__class__, self).__init__(workdir,
                                             magnetizingCurves, magnets,
                                             result_func,
                                             template_name='sampling')


class Moga(Dakota):
    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, result_func=None):
        super(self.__class__, self).__init__(workdir,
                                             magnetizingCurves, magnets,
                                             result_func,
                                             template_name='moga')


class PsuadeMoat(Dakota):
    def __init__(self, workdir,
                 magnetizingCurves=None, magnets=None, result_func=None):
        super(self.__class__, self).__init__(workdir,
                                             magnetizingCurves, magnets,
                                             result_func,
                                             template_name='psuade_moat')
