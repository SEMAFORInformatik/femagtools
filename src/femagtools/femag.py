"""running FEMAG

"""
import subprocess
import pathlib
import os
import glob
import logging
try:
    import zmq
except ImportError:
    pass
import sys
import json
import io
import time
import re
import threading
import femagtools.model
import femagtools.magnet
import femagtools.conductor
import femagtools.windings
import femagtools.mcv
import femagtools.asm
import femagtools.airgap as ag
import femagtools.fsl
import femagtools.config
import femagtools.ecloss
import femagtools.forcedens
import femagtools.zmq

from femagtools import ntib


logger = logging.getLogger(__name__)


class FemagError(Exception):
    pass


def handle_process_output(filedes, outfile, log):
    """read from file descriptor and direct lines to logger and outfile"""
    with open(outfile, 'wb') as fp:
        for line in filedes:
            fp.write(line)
            if log:
                if (b'' == line or
                    b'\x1b' in line or  # ignore terminal escape seq
                    b'\n' == line or
                    b'Start:' in line or
                    b'Stop:' in line or
                        b'Elapsed time:' in line):
                    continue
                if line:
                    try:
                        logger.info(" > %s",
                                    line.decode().strip())
                    except UnicodeDecodeError:
                        # ignore
                        pass


def set_magnet_properties(model, simulation, magnets):
    """set temperature adapted magnet properties"""
    if not hasattr(model, 'magnet'):
        return
    if isinstance(simulation, dict):
        if 'hc_min' in simulation:  # backward compatibility
            model.hc_min = simulation['hc_min']
        if 'magn_temp' not in simulation:
            return
        magn_temp = simulation['magn_temp']
    else:
        if hasattr(simulation, 'hc_min'):  # backward compatibility
            model.hc_min = simulation.hc_min
        if not hasattr(simulation, 'magn_temp'):
            return
        magn_temp = simulation.magn_temp
    if 'material' not in model.magnet:
        return
    try:
        material = model.magnet['material']
        magnetMat = magnets.find(material)
        if magnetMat:
            model.magnet['temp_prop'] = {
                'remanenc': magnetMat.get('remanenc', 1.2),
                'magntemp': magn_temp
            }
            relperm = magnetMat.get('relperm', 1.05)
            tempcoefmuer = 1
            tempcoefbr = magnetMat.get('temcoefbr', 0)
            tempcoefhc = magnetMat.get('temcoefhc', 0)
            if tempcoefbr and tempcoefhc:
                tempcoefmuer = (tempcoefbr-tempcoefhc)/(1+tempcoefhc*(magn_temp-20))
            if 'temcoefmuer' in magnetMat:
                tempcoefmuer = magnetMat['temcoefmuer']
            hcj = magnetMat.get('HcJ', 0)
            if hcj:
                model.hc_min = -hcj * (1+tempcoefhc*(magn_temp-20))
            model.magnet['temp_prop']['relperm'] = \
                (1+tempcoefmuer*(magn_temp-20))*relperm
            if tempcoefbr:
                model.magnet['temp_prop']['temcoefbr'] = tempcoefbr
            if tempcoefhc:
                model.magnet['temp_prop']['temcoefhc'] = tempcoefhc

    except AttributeError:
        pass

class BaseFemag(object):
    def __init__(self, workdir, cmd, magnetizingCurves, magnets, condMat,
                 templatedirs=[]):
        self.workdir = workdir
        self.magnets = []
        self.magnetizingCurves = []
        self.condMat = []
        self.templatedirs = templatedirs
        if cmd:
            self.cmd = cmd
        else:
            self.cmd = ''

        if magnetizingCurves:
            if isinstance(magnetizingCurves,
                          femagtools.mcv.MagnetizingCurve):
                self.magnetizingCurves = magnetizingCurves
            else:
                self.magnetizingCurves = femagtools.mcv.MagnetizingCurve(
                    magnetizingCurves)

        if magnets:
            if isinstance(magnets, femagtools.magnet.Magnet):
                self.magnets = magnets
            else:
                self.magnets = femagtools.magnet.Magnet(magnets)

        if condMat:
            if isinstance(condMat, femagtools.conductor.Conductor):
                self.condMat = condMat
            else:
                self.condMat = femagtools.conductor.Conductor(condMat)

    def copy_winding_file(self, name, wdg):
        wdg.write(name, self.workdir)

    def copy_magnetizing_curves(self, model, dir=None, recsin='', feloss=''):
        """extract mc names from model and write files into workdir or dir if given

        Return:
            list of extracted mc names (:obj:`list`)
        """
        dest = dir if dir else self.workdir
        if isinstance(feloss, (int, float)):
            try:
                feloss = {1: 'jordan', 11: 'bertotti'}[int(feloss)]
            except KeyError:
                feloss = ''
        return [self.magnetizingCurves.writefile(m[0], dest,
                                                 fillfac=m[1],
                                                 recsin=recsin,
                                                 feloss=feloss)
                for m in model.set_magcurves(
            self.magnetizingCurves, self.magnets)]

    def create_wdg_def(self, model):
        name = 'winding'
        w = femagtools.windings.Winding(
            dict(
                Q=model.stator['num_slots'],
                p=model.poles//2,
                m=len(model.winding['wdgdef']),
                windings=model.winding['wdgdef']))
        self.copy_winding_file(name, w)
        return name

    def create_fsl(self, machine, simulation):
        """create list of fsl commands
        Args:
          machine: dict of fe machine model (stator, rotor/magnet, winding)
          simulation: dict of FE simulation (can be empty)
        """
        self.model = femagtools.model.MachineModel(machine)
        self.modelname = self.model.name
        recsin = ''
        feloss = ''
        if simulation:
            recsin = simulation.get('recsin', '')
            feloss = simulation.get('calc_fe_loss', '')
        self.copy_magnetizing_curves(self.model, recsin=recsin, feloss=feloss)
        try:
            if 'wdgdef' in self.model.winding:
                self.model.winding['wdgfile'] = self.create_wdg_def(
                    self.model)
        except AttributeError:
            pass
        builder = femagtools.fsl.Builder(self.templatedirs)
        if simulation:
            if 'num_par_wdgs' not in simulation:
                try:
                    num_par_wdgs = self.model.winding['num_par_wdgs']
                    simulation['num_par_wdgs'] = num_par_wdgs
                except:
                    pass
            set_magnet_properties(self.model, simulation, self.magnets)
            return builder.create(self.model, simulation,
                                  self.magnets, self.condMat)
        return builder.create_model(self.model,
                                    self.magnets,
                                    self.condMat) + ['save_model("cont")']

    def get_log_value(self, pattern, modelname='FEMAG-FSL.log'):
        result = []
        pat = re.compile(pattern)
        with open(os.path.join(self.workdir, 'FEMAG-FSL.log')) as f:
            for line in f:
                m = pat.search(line)
                if m:
                    result.append(float(line.split(':')[-1].split()[0]))
        return result

    def get_bch_file(self, modelname, offset=0):
        return self.get_result_file(modelname, 'B(AT)?CH', offset)

    def get_plt_file(self, modelname, offset=0):
        return self.get_result_file(modelname, 'PLT0', offset)

    def get_asm_file(self, modelname, offset=0):
        return self.get_result_file(modelname, 'ASM', offset)

    def get_ts_files(self, modelname, offset=0):
        return self.get_result_file(modelname, 'TS', offset)

    def get_result_file(self, modelname, ext, offset=0):
        """return latest result (bch, asm) file (if any)"""
        filelist = glob.glob(os.path.join(
            self.workdir, modelname+'_[0-9][0-9][0-9].*'))
        if (filelist):
            return sorted([f for f in filelist
                           if re.search(ext, f.split('.')[-1])])[-1-offset]
        return ''

    def get_result_file_list(self, modelname, ext):
        """return list of result (me) file (if any)"""
        filelist = sorted(glob.glob(os.path.join(
            self.workdir, modelname+'_[0-9][0-9][0-9].'+ext)))
        if len(filelist) > 0:
            return filelist
        return []

    def read_asm(self, modelname=None, offset=0):
        "read most recent ASM file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()

        asmfile = self.get_asm_file(modelname, offset)
        if asmfile:
            logger.info("Read ASM {}".format(asmfile))
            return femagtools.asm.read(asmfile)
        return {}

    def read_ts(self, modelname=None):
        "read most recent TS files and return result"
        import femagtools.ts
        if not modelname:
            modelname = self._get_modelname_from_log()

        logger.info("Read TS {}".format(modelname))
        return femagtools.ts.read_st(self.workdir, modelname)

    def read_modal(self, modelname=None):
        "read modal analysis output and return result"
        import femagtools.me
        filelist = self.get_result_file_list(modelname + '_Mode', 'txt')
        psfiles = self.get_result_file_list(modelname + '_Mode', 'ps')
        if len(filelist) > 0:
            logger.info("Read Eigenvectors {}".format(modelname))
            return femagtools.me.get_eigenvectors(filelist, psfiles)
        else:
            return ''

    def read_bch(self, modelname=None, offset=0):
        "read most recent BCH/BATCH file and return result"
        # read latest bch file if any
        if not modelname:
            modelname = self._get_modelname_from_log()

        result = femagtools.bch.Reader()
        bchfile = self.get_bch_file(modelname, offset)
        if bchfile:
            logger.info("Read BCH {}".format(bchfile))
            with io.open(bchfile, encoding='latin1',
                         errors='ignore') as f:
                result.read(f)
        return result

    def read_forcedens(self, modelname=None, offset=0):
        "read most recent PLT0 file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()

        pltfile = self.get_plt_file(modelname, offset)
        if pltfile:
            logger.info("Read PLT0 {}".format(pltfile))
            return femagtools.forcedens.read(pltfile)
        return None

    def read_isa(self, modelname=None):
        "read most recent I7/ISA7 file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()
        import femagtools.isa7
        return femagtools.isa7.read(os.path.join(self.workdir, modelname))

    def read_nc(self, modelname=None):
        "read most recent NC file and return result"
        if not modelname:
            modelname = self._get_modelname_from_log()
        import femagtools.nc
        return femagtools.nc.read(os.path.join(self.workdir, modelname))

    def read_los(self, modelname=None):
        "read most recent LOS file and return result"
        # read latest los file if any
        if not modelname:
            modelname = self._get_modelname_from_log()

        losfile_list = sorted(glob.glob(os.path.join(
            self.workdir, modelname+'_[0-9][0-9][0-9].LOS')))
        if len(losfile_list) > 0:
            return ntib.read_los(losfile_list[-1])

        return dict()

    def read_hsn(self, modelname=None):
        import numpy as np
        "read heat network result"
        _map = {
            "StZa": "plfe1",
            "outs": "-",
            "StJo": "plfe1",
            "Slot": "-",
            "Shaf": "-",
            "Iron": "plfe2",
            "PMag": "plmag",
            "PMag_1": "plmag",
            "PMag_2": "plmag",
            "PMag_3": "plmag",
            "PMag_4": "plmag",
            "W1  ":   "plcu1",
            "W2  ":   "plcu1",
            "W3  ":   "plcu1"
        }
        if not modelname:
            modelname = self._get_modelname_from_log()
        hsn_list = sorted(glob.glob(os.path.join(
            self.workdir, modelname+'.hsn')))
        with open(hsn_list[-1], 'r') as f:
            hsn_data = json.load(f)

        # area calculation
        nc_file = self.read_nc(modelname)
        slot_area = 0.0
        wdg_area = []
        for i in nc_file.windings:
            for j in i.subregions:
                wdg_area.append(j.area())

        slot_area = np.sum(wdg_area).item()/3
        magnet_area = 0.0
        num_sreg_mag = 0
        area = dict()
        for i in nc_file.subregions:
            if i.name in ("StZa", "StJo", "Iron"):
                area[i.name] = i.area().item()
            elif i.name == "PMag":
                magnet_area += i.area().item()
                num_sreg_mag += 1
            else:
                pass

        area['PMag'] = magnet_area/num_sreg_mag
        for i in ('W1', 'W2', 'W3'):
            area[i] = slot_area

        pmag_index = []
        if "Nodes" in hsn_data:
            for k ,i in enumerate(hsn_data['Nodes']):
                i.update({"mass": i['weight'], "losses": _map[i['Name']]})
                if "PMag" in i['Name']:
                    pmag_index.append(k)
                if i['Name'].strip() in area.keys():
                    i.update({"area": area[i['Name'].strip()]})
            if pmag_index:
                for i in range(len(pmag_index)):
                    hsn_data["Nodes"][pmag_index[i]]['Name'] = f"PMag_{i+1}"
        with open(hsn_list[-1], 'w') as f:
            json.dump(hsn_data, f)
        return nc_file

    def _get_modelname_from_log(self):
        """
        Read the modelname from the Femag Log file
        """
        try:
            with open(os.path.join(self.workdir, 'FEMAG-FSL.log')) as f:
                for l in f:
                    if l.startswith('New model') or l.startswith('Load model'):
                        return l.split('"')[1]
        except FileNotFoundError:
            pass

        return list(pathlib.Path(self.workdir).glob('*.PROT'))[0].stem

    def readResult(self, machine, simulation, bch=None):
        if simulation:
            if simulation['calculationMode'] == "fieldcalc":
                nc = self.read_nc()
                pmod = nc.poles_sim
                r = {'airgap': ag.read(
                    os.path.join(self.workdir, 'bag.dat'), pmod=pmod)}
                if 'plots' in simulation:
                    if 'field_lines' in simulation['plots']:
                        r['field_lines'] = os.path.join(
                            self.workdir, 'field.svg')
                return r
            if simulation['calculationMode'] == "pm_sym_loss":
                return self.read_los(self.modelname)

            if simulation['calculationMode'] == 'asyn_motor':
                return self.read_asm(self.modelname)

            if simulation['calculationMode'] == 'calc_field_ts':
                return self.read_ts(self.modelname)

            if simulation['calculationMode'] == 'modal_analysis':
                return self.read_modal(self.modelname)

            if simulation['calculationMode'] == 'therm-dynamic':
                temp = [[float(n) for n in l.split()]
                        for l in (pathlib.Path(self.workdir) /
                                  'temperature.dat').read_text().split('\n') if l]
                ttemp = list(zip(*temp))
                return {'t': ttemp[0], 'temperature': ttemp[1]}

            if simulation['calculationMode'] == 'hsn':
                model = None
                try:
                    model = self.read_hsn()
                except:
                    pass
                if model is None:
                    model = self.read_nc()
                return model.get_minmax_temp()

            if not bch:
                bch = self.read_bch(self.modelname)
            if simulation['calculationMode'] == 'pm_sym_fast' or \
                simulation['calculationMode'] == 'torq_calc':
                if simulation.get('shortCircuit', False):
                    from .shortcircuit import shortcircuit
                    set_magnet_properties(self.model, simulation, self.magnets)
                    bch.scData = shortcircuit(self, machine, bch, simulation)
                    #bch.torque += bchsc.torque
                    #bch.demag += bchsc.demag

            if 'airgap_induc' in simulation:
                try:
                    pmod = bch.machine['p_sim']
                except KeyError:
                    pmod = 0
                bch.airgap = ag.read(os.path.join(self.workdir, 'bag.dat'),
                                     pmod=pmod)

            if simulation.get('magnet_loss', False):
                logger.info('Evaluating magnet losses...')
                ops = range(len(bch.torque))
                ncf = pathlib.Path(self.workdir) / self.modelname
                m = femagtools.ecloss.MagnLoss(
                    nc=femagtools.nc.read(ncf), ibeta=ops)
                try:
                    # change from ialh to ialh2: since v1.8.1
                    magn_losses = m.calc_losses_ialh2()
                except:
                    magn_losses = list(range(len(ops)))

                if len(ops) != len(bch.losses):
                    magn_losses.insert(0, magn_losses[0])
                try:
                    for i in range(len(bch.losses)):
                        bch.losses[i].update({"magnetH": magn_losses[i]})
                except:
                    pass
                # pass losses to bch object for th usage
                try:
                    bch.magnet_loss_th = m.th_loss
                except:
                    pass
            try:
                if hasattr(self, 'dy2'):
                    setattr(bch, 'dy2', self.dy2)
            except:
                pass
            return bch


class Femag(BaseFemag):
    """Invoke and control execution of FEMAG

    Args:
        workdir: name of working directory
        cmd: name of femag program (default is config.executable)
        magnetizingCurves: collection of lamination material curves or name of directory
        magnets: collection of magnet material
        condMat: collection of conductor material
        templatedirs: (list) names of directories that include mako files as fsl templates

    """

    def __init__(self, workdir, cmd=None, templatedirs=[],
                 magnetizingCurves='.', magnets=None, condMat=[]):
        super(self.__class__, self).__init__(workdir, cmd,
                                             magnetizingCurves, magnets, condMat,
                                             templatedirs=templatedirs)

    def run(self, filename, options=['-b'], fsl_args=[],
            stateofproblem='mag_static'):
        """invoke FEMAG in current workdir

        Args:
            filename: name of FSL file to execute
            options: list of FEMAG options
            fsl_args: list of FSL argument options
            stateofproblem: (str) one of config.state_of_problem_set
        Raises:
            FemagError
        """
        if self.cmd:
            cmd = self.cmd
        else:
            cmd = femagtools.config.get_executable(stateofproblem)
        if (cmd.find('wfemag') > -1 and
            '-b' in options and
                '-m' not in options):
            options.insert(0, '-m')
        args = [cmd] + options + [filename] + fsl_args

        basename = pathlib.Path(filename).name
        outname = os.path.join(self.workdir, basename+'.out')
        errname = os.path.join(self.workdir, basename+'.err')
        with open(errname, 'w') as err:
            logger.info('invoking %s', ' '.join([str(a) for a in args]))
            proc = subprocess.Popen(
                args,
                stdout=subprocess.PIPE, stderr=err, cwd=self.workdir)
            stdoutthr = threading.Thread(target=handle_process_output,
                                         args=(proc.stdout, outname, True))
            stdoutthr.start()
        proc.wait()
        stdoutthr.join()
        errs = []
        # print femag output
        with io.open(outname, encoding='latin1', errors='ignore') as outfile:
            errLine = False
            for l in outfile:
                if l.find('ERROR') > -1:
                    errs.append(l.strip())
                    errLine = True
                elif errLine and l.startswith(' '):  # additional error line
                    errs.append(l)
                else:
                    errLine = False

        rc = proc.returncode
        logger.info("%s exited with returncode %d (num errs=%d)",
                    cmd, rc, len(errs))
        if rc != 0 or errs:
            with io.open(errname, encoding='latin1',
                         errors='ignore') as errfile:
                for l in errfile:
                    errs.append(l.strip())
            errs.insert(0, 'Exit code {}'.format(rc))
            raise FemagError("\n".join(errs))

    def cleanup(self):
        "removes all created files in workdir"
        if not os.path.exists(self.workdir):
            return
        cleanfiles = ('*.B*CH', '*.I*7-*', '*.A*7-*', '*.nc-*',
                      '*.dat', '*.LOS', '*.svg', '*.png', '*.hxy')
        # '*.TMC','*.TMO', '*.PROT'):
        for p in cleanfiles:
            for f in glob.glob(os.path.join(self.workdir, p)):
                os.remove(f)

    def __call__(self, machine, simulation={},
                 options=['-b'], fsl_args=[]):
        """setup fsl file, run calculation and return
        BCH, ASM, TS or LOS results if any."""
        fslfile = 'femag.fsl'
        with open(os.path.join(self.workdir, fslfile), 'w') as f:
            f.write('\n'.join(self.create_fsl(machine,
                                              simulation)))
        if simulation:
            stateofproblem = simulation.get('stateofproblem', 'mag_static')
            if 'eccentricity' in simulation:
                #assert (self.model.stator['num_slots']
                #        == self.model.stator['num_slots_gen'])
                for k in ('bore_diam', 'airgap'):
                    if k not in simulation['eccentricity']:
                        simulation['eccentricity'][k] = self.model.get(k)

            if 'poc' in simulation:
                with open(os.path.join(self.workdir,
                                       simulation['pocfilename']), 'w') as f:
                    f.write('\n'.join(simulation['poc'].content()))
            if simulation['calculationMode'] == "pm_sym_loss":
                with open(os.path.join(self.workdir,
                                       self.modelname+'.ntib'), 'w') as f:
                    f.write('\n'.join(ntib.create(
                        simulation['speed'],
                        simulation['current'],
                        simulation['angl_i_up'])))
                # TODO: add r1, m
            if simulation['calculationMode'] == 'therm-dynamic':
                load = simulation['load']
                (pathlib.Path(self.workdir)/'load.csv').write_text(
                    '\n'.join([','.join([f"{load[k][i]}"
                                         for k in ('t', 'n', 'T', 'i1', 'beta')])
                               for i in range(0, len(load['t']))]))
        else:
            stateofproblem = 'mag_static'

        self.run(fslfile, options, fsl_args, stateofproblem=stateofproblem)

        try:
            setattr(self, "dy2", machine['stator']['dy2'])
        except:
            pass

        if simulation:
            return self.readResult(machine, simulation)

        return {'status': 'ok', 'message': self.modelname,
                'model': self.model.props()}

class FemagTask(threading.Thread):
    def __init__(self, port, args, workdir, logdir):
        threading.Thread.__init__(self)
        self.port = port
        self.args = args + [str(self.port)]
        self.proc = None
        self.returncode = None
        self.workdir = workdir
        self.logdir = logdir

    def run(self):
        logger.info("femag is ready on port %d workdir %s",
                    self.port, self.workdir)
        outname = os.path.join(self.logdir, f'femag-{self.port}.out')
        errname = os.path.join(self.logdir, f'femag-{self.port}.err')
        with open(outname, 'w') as out, open(errname, 'w') as err:
            self.proc = subprocess.Popen(
                self.args,
                stdout=out, stderr=err, cwd=self.workdir)

        self.returncode = self.proc.wait()


class ZmqFemag(BaseFemag):
    """Invoke and control execution of FEMAG with ZeroMQ

    Args:
        port: port number of req socket
        host: hostname (ip addr)
        workdir: name of working directory
        logdir: name of logging directory (default is workdir/log)
        cmd: name of femag program
        templatedirs: (list) names of directories that include mako files as fsl templates
    """

    def __init__(self, port, host='localhost', workdir='', logdir='',
                 cmd=None,
                 magnetizingCurves=None, magnets=None, condMat=[],
                 templatedirs=[]):
        super(self.__class__, self).__init__(
            workdir, cmd,
            magnetizingCurves, magnets, condMat,
            templatedirs=templatedirs)
        self.host = host
        self.port = port
        self.femaghost = ''
        self.logdir = logdir if logdir else os.path.join(workdir, 'log')
        if not os.path.exists(self.logdir):
            os.makedirs(self.logdir)
        self.request_socket = None
        self.subscriber = None

    def close(self):
        if self.subscriber:
            self.subscriber.stop()

        if self.request_socket:
            logger.debug("close request_socket")
            self.request_socket.close()
            self.request_socket = None

        logger.debug("done")

    def __del__(self):
        try:
            self.close()
        except AttributeError:
            pass
        logger.debug("Destructor ZmqFemag")

    def stopFemag(self):
        if self.femagTask and self.femagTask.proc:
            logger.info("stopFemagTask")
            self.femagTask.proc.kill()
            self.femagTask.proc.wait()
            self.femagTask = None
            self.request_socket = None
            self.subscriber.stop()
            self.subscriber = None
            logger.info("stopFemagTask Done")
        else:
            logger.warning("stopFemag not implemented")

    def __req_socket(self):
        """returns a new request client"""
        context = zmq.Context.instance()
        request_socket = context.socket(zmq.REQ)
        url = 'tcp://{0}:{1}'.format(
            self.host, self.port)
        logger.debug("connect req socket %s", url)
        request_socket.connect(url)
        return request_socket

    def subscribe(self, notify):
        """attaches a notify function"""
        logger.info("Subscribe on '%s' port %d", self.femaghost, self.port+1)
        femagtools.zmq.SubscriberTask.clear()
        if self.subscriber is None:
            # progress/xyplot at a configured timestep published
            header = [b'progress', b'xyplot', b'license']
            self.subscriber = femagtools.zmq.SubscriberTask(
                port=self.port+1, host=self.femaghost, notify=notify, header=header)
            self.subscriber.start()
        else:
            # reattach?
            self.subscriber.notify = notify
        return

    def __is_running(self):
        try:
            return self.femagTask.proc and self.femagTask.returncode is None
        except:
            pass
        return False

    def send_request(self, msg, pub_consumer=None, timeout=None):
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        if not msg:
            return [b'{"status":"ignored",{}}']
        if timeout:
            self.request_socket.setsockopt(zmq.RCVTIMEO, timeout)
            self.request_socket.setsockopt(zmq.LINGER, 0)
        else:
            self.request_socket.setsockopt(zmq.RCVTIMEO, -1)  # default value

        for m in msg[:-1]:
            self.request_socket.send_string(m, flags=zmq.SNDMORE)
        if isinstance(msg[-1], list):
            self.request_socket.send_string('\n'.join(msg[-1]))
        else:
            self.request_socket.send_string(msg[-1])

        errmsg = ''
        max_num_again = 2
        again = 0
        while True:
            try:
                return self.request_socket.recv_multipart()
            except zmq.error.Again as e:
                again += 1
                # logger.exception("send_request")
                errmsg = str(e)
                logger.warning("send_request: %s Message %s, host: %s port: %s", str(
                    e), msg, self.host, self.port)
                if again >= max_num_again:
                    break
                continue
        logger.info("send_request failed %s", errmsg.encode())
        return [b'{"status":"error", "message":"' + errmsg.encode() + b'"}']

    def send_fsl(self, fsl, timeout=None):
        """sends FSL commands in ZMQ mode and blocks until commands are processed

        Args:
            fsl: string (or list) of FSL commands
            timeout: The maximum time (in milliseconds) to wait for a response

        Return:
            status
        """
        header = 'FSL'
        self.request_socket.send_string(header, flags=zmq.SNDMORE)

        if isinstance(fsl, list):
            self.request_socket.send_string('\n'.join(fsl))
        else:
            self.request_socket.send_string(fsl)
        logger.debug("Sent fsl wait for response")

        if timeout:
            self.request_socket.setsockopt(zmq.RCVTIMEO, timeout)
            self.request_socket.setsockopt(zmq.LINGER, 0)
        else:
            self.request_socket.setsockopt(zmq.RCVTIMEO, -1)  # default value
        import datetime
        startTime = datetime.datetime.now() if timeout else None
        while True:
            try:
                response = self.request_socket.recv_multipart()
                # NOTE: femag encoding is old school
                return [s.decode('latin1') for s in response]
            except zmq.error.Again as e:
                logger.info("Again [%s], timeout: %d", str(e), timeout)
                if not startTime:
                    continue

                diffTime = datetime.datetime.now() - startTime
                logger.debug("Diff msec[%d]", diffTime.microseconds)
                if diffTime.microseconds < timeout:
                    continue

                logger.info("ALERT not running close socket")
                self.close()
                return ['{"status":"error", "message":"Femag is not running"}', '{}']

            except Exception as e:
                logger.exception("send_fsl")
                logger.error("send_fsl: %s", str(e))
                if timeout:  # only first call raises zmq.error.Again
                    return ['{"status":"error", "message":"Femag is not running"}', '{}']
                msg = json.dumps(str(e))
                return ['{"status":"error", "message":'+msg+'}', '{}']

    def run(self, options=['-b'], restart=False, procId=None,
            stateofproblem='mag_static'):  # noqa: C901
        """invokes FEMAG in current workdir and returns pid

        Args:
            options: list of FEMAG options
            stateofproblem: str one of config.executable
        """
        if self.__is_running():
            if restart:
                logger.info("must restart")
                self.quit(True)
                self.femagTask.join()
                logger.info("Stopped procId: %s", self.femagTask.proc.pid)
            else:
                return self.femagTask.proc.pid

        if self.cmd:
            cmd = self.cmd
        else:
            cmd = femagtools.config.get_executable(stateofproblem)
        if (cmd.find('wfemag') > -1 and
            '-b' in options and
                '-m' not in options):
            options.insert(0, '-m')

        args = [cmd] + options
        self.femagTask = FemagTask(self.port, args, self.workdir, self.logdir)
        self.femagTask.start()
        if not self.request_socket:
            self.request_socket = self.__req_socket()

        # check if mq is ready for listening
        lcount = 300
        for t in range(lcount):
            time.sleep(0.1)
            if self.__is_running():
                if self.femagTask.proc:
                    logger.info("femag (pid: '{}') is listening".format(
                        self.femagTask.proc.pid))
                break

        return self.femagTask.proc.pid if self.femagTask.proc else 0

    def quit(self, save_model=False):
        """terminates femag"""

        if not self.__is_running():
            logger.info("Femag already stopped")
            return

        if self.subscriber:
            self.subscriber.stop()
            self.subscriber = None

        # send exit flags
        f = '\n'.join(['exit_on_end = true',
                       'exit_on_error = true'])
        response = self.send_fsl(f)

        # send quit command
        try:
            response = [r.decode('latin1')
                        for r in self.send_request(
                ['CONTROL', 'quit'], timeout=2000)]
        except Exception as e:
            logger.error("Femag Quit zmq message %s", e)

        logger.debug("Sent QUIT to femag %s", response)
        # if query, send a answer
        obj = json.loads(response[0])
        logger.debug("status: {}".format(obj['status']))
        if obj['status'] == 'Query':
            logger.info('query: %s => %s',
                        obj['message'], 'saved' if save_model else 'not saved')

            # Only send one msg
            response = self.request_socket.send_string(
                'Ok' if save_model else 'Cancel')

    def upload(self, files):
        """upload file or files
        returns number of transferred bytes
        (FEMAG 8.5 Rev 3282 or greater only)
        """
        ret = []
        total = 0
        chunk_size = 20*1024
        fnames = files
        if isinstance(files, str):
            fnames = [files]
        for fn in fnames:
            self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
            basename = os.path.basename(fn)
            logger.info("upload %s --> %s", fn, basename)
            self.request_socket.send_string('upload = {}'.format(basename),
                                            flags=zmq.SNDMORE)
            with open(fn, mode="rb") as file:
                while True:
                    data = file.read(chunk_size)
                    if not data:
                        break
                    more = 0 if len(data) < chunk_size else zmq.SNDMORE
                    self.request_socket.send(data, flags=more)
                    total += len(data)
            ret += [s.decode('latin1')
                    for s in self.request_socket.recv_multipart()]
        return ret

    def copyfile(self, filename, dirname):
        """copy filename to dirname (FEMAG 9.2)"""
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
        self.request_socket.send_string(f'copyfile {filename} {dirname}')
        return [r.decode() for r in self.request_socket.recv_multipart()]

    def change_case(self, dirname):
        """change case to dirname (FEMAG 9.2)"""
        logger.info("change_case to :  %s", dirname)
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
        self.request_socket.send_string(f'casedir {dirname}')
        return [r.decode() for r in self.request_socket.recv_multipart()]

    def delete_case(self, dirname):
        """delete case dir (FEMAG 9.2)"""
        if not self.request_socket:
            self.request_socket = self.__req_socket()
        self.request_socket.send_string('CONTROL', flags=zmq.SNDMORE)
        self.request_socket.send_string(f'casedir -{dirname}')
        return [r.decode() for r in self.request_socket.recv_multipart()]

    def clear(self, timeout=2000):
        """clear lua script session"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'clear'],
                                           timeout=timeout)]

    def cleanup(self, timeout=2000):
        """remove all FEMAG files in working directory
        (FEMAG 8.5 Rev 3282 or greater only)"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'cleanup'],
                                           timeout=timeout)]

    def release(self):
        """signal finish calculation task to load balancer to free resources
        (Docker Cloud environment only)
        """
        return [r.decode('latin1')
                for r in self.send_request(['close'], timeout=1000)]

    def info(self, timeout=2000):
        """get various resource information
        (FEMAG 8.5 Rev 3282 or greater only)"""
        response = [r.decode('latin1')
                    for r in self.send_request(['CONTROL', 'info'],
                                               timeout=timeout)]
        status = json.loads(response[0])['status']
        if status == 'ok':
            try:
                self.femaghost = json.loads(response[1])['addr']
                logger.info("Set femaghost %s", self.femaghost)
            except KeyError:
                # ignore femaghost setting if no addr is returned (Windows only)
                pass
        return response

    def publishLevel(self, level):
        """set publish level"""
        return [r.decode('latin1')
                for r in self.send_request(['CONTROL', 'publish = {}'
                                            .format(level)], timeout=10000)]

    def getfile(self, filename=''):
        """get file (FEMAG 8.5 Rev 3282 or greater only)"""
        response = self.send_request(
            ['CONTROL', f'getfile = {filename}'], timeout=1000)
        return [response[0].decode('latin1'),
                response[1] if len(response) > 1 else b'']

    def exportsvg(self, fslcmds, timeout=10000):
        """get svg format from fsl commands (if any graphic created)
        (since FEMAG 8.5 Rev 3343) """
        response = self.send_request(['SVG', fslcmds], timeout=timeout)
        try:
            rc = json.loads(response[0].decode('latin1'))
            if rc['status'] == 'ok':
                return self.getfile(rc['result_file'][0])
        except json.decoder.JSONDecodeError:
            logger.warning(response[0])
        return [s.decode('latin1') for s in response]

    def exportmesh(self, fslcmds, timeout=120000):
        """get svg format (with mesh) from fsl commands (if any graphic created)
        (since FEMAG v2024.1-17) """
        response = self.send_request(['SVG+Mesh', fslcmds], timeout=timeout)
        try:
            rc = json.loads(response[0].decode('latin1'))
            if rc['status'] == 'ok':
                return self.getfile(rc['result_file'][0])
        except json.decoder.JSONDecodeError:
            logger.warning(response[0])
        return [s.decode('latin1') for s in response]

    def airgap_flux_density(self, pmod):
        # try to read bag.dat
        agr = self.getfile("bag.dat")
        status = json.loads(agr[0])['status']
        if status == 'ok':
            datfile = os.path.join(self.workdir, 'bag.dat')
            with open(datfile, 'wb') as bagfile:
                bagfile.write(agr[1])
            agi = ag.read(datfile, pmod)
        else:
            import numpy as np
            # try to read model file (TODO download with getfile)
            nc = self.read_nc()
            ag_elmnts = nc.airgap_center_elements
            logger.info("Airgap elements %d scale_factor %f",
                        len(ag_elmnts), nc.scale_factor())
            if len(ag_elmnts) < 1:
                raise ValueError("Missing airgap elements")
            scf = 360/nc.scale_factor()/ag_elmnts[-1].center[0]
            pos = np.array([e.center[0]*scf for e in ag_elmnts])
            bxy = np.array([e.flux_density() for e in ag_elmnts]).T
            if np.max(bxy[0]) > np.max(bxy[1]):
                agi = ag.fft(pos, bxy[0], pmod)
            else:
                agi = ag.fft(pos, bxy[1], pmod)
        return dict(Bamp=agi['Bamp'],
                    phi0=agi['phi0'],
                    angle=agi['pos'],
                    angle_fft=agi['pos'],
                    B=agi['B'],
                    B_fft=agi['B_fft'])

    def interrupt(self):
        """send push message to control port to stop current calculation"""
        context = zmq.Context.instance()
        if not self.femaghost:
            self.femaghost = '127.0.0.1'
        ctrl = context.socket(zmq.PUSH)
        ctrl.connect('tcp://{0}:{1}'.format(
            self.femaghost, self.port+2))
        logger.info("Interrupt %s", self.femaghost)
        ctrl.send_string('interrupt')
        ctrl.close()
        femagtools.zmq.SubscriberTask.clear()

    def copy_winding_file(self, name, wdg):
        wdg.write(name, self.workdir)
        self.upload(os.path.join(self.workdir, name+'.WID'))

    def copy_magnetizing_curves(self, model, dir=None, recsin='', feloss=''):
        """extract mc names from model and write files into workdir or dir if given
           and upload to Femag

        Return:
            list of extracted mc names (:obj:`list`)
        """
        dest = dir if dir else self.workdir
        mc_names = [m for m in model.set_magcurves(
            self.magnetizingCurves, self.magnets)]
        if isinstance(feloss, (int, float)):
            try:
                feloss = {1: 'jordan', 11: 'bertotti'}[int(feloss)]
            except KeyError:
                feloss = ''
        for m in mc_names:
            f = self.magnetizingCurves.writefile(m[0], dest,
                                                 fillfac=m[1],
                                                 recsin=recsin,
                                                 feloss=feloss)
            self.upload(os.path.join(dest, f))
        return mc_names

    def __call__(self, machine, simulation):
        """setup fsl file, run calculation and return BCH results
        Args:
          machine: dict with machine parameters or name of model
          simulation; dict with simulation parameters

        Raises:
           FemagError
        """
        if isinstance(machine, str):
            modelpars = dict(name=machine)
        else:
            modelpars = machine
        if 'exit_on_end' not in modelpars:
            modelpars['exit_on_end'] = 'false'
        if 'exit_on_error' not in modelpars:
            modelpars['exit_on_error'] = 'false'
        response = self.send_fsl(['save_model("close")'] +
                                 self.create_fsl(modelpars,
                                                 simulation))
        r = json.loads(response[0])
        if r['status'] != 'ok':
            raise FemagError(r['message'])

        if not simulation:
            model = self.model.props()
            model['name'] = machine['name']
            return [r, model]

        result_file = r['result_file'][0]
        if simulation['calculationMode'] == "pm_sym_loss":
            return self.read_los(self.modelname)

        if simulation['calculationMode'] == "asyn_motor":
            results = self.read_asm(self.modelname)
            if 'airgap_induc' in simulation:
                try:
                    pmod = results['p_gen']
                except KeyError:
                    pmod = 0
                bagdat = os.path.join(self.workdir, 'bag.dat')
                results['airgap'] = ag.read(bagdat, pmod=pmod)
            return results

        if simulation['calculationMode'] == "calc_field_ts":
            return self.read_ts(self.modelname)

        status, content = self.getfile(result_file)
        r = json.loads(status)
        if r['status'] == 'ok':
            bch = femagtools.bch.Reader()
            bch.read(content.decode('latin1'))
            bch = self.readResult(simulation, bch)
            return bch
        raise FemagError(r['message'])
