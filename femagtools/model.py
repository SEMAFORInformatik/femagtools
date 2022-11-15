# -*- coding: utf-8 -*-
"""
    femagtools.model
    ~~~~~~~~~~~~~~~~

    Managing model parameters



"""
import logging
import string
import numpy as np

logger = logging.getLogger(__name__)
#
# Set of legal model name chars
# 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ
MODELNAME_CHARS = string.printable[:63] + "-_äöüéè"
# maximum name length
MAX_MODEL_NAME_LEN = 255


def movesteps(nodes):
    """returns list of move steps

    Args:
       nodes: number of nodes in airgap
    """
    return [nodes // w + 1 for w in range(1, nodes//2) if nodes % w == 0]


class MCerror(Exception):
    pass


class Model(object):
    def __init__(self, parameters):
        if isinstance(parameters, dict):
            for k in parameters.keys():
                setattr(self, k, parameters[k])

    def set_value(self, name, value, p=None):
        """set value of parameter identified by name

        Args:
            name: name of parameter
            value: value to be assigned to parameter
        """
        if isinstance(name, str):
            setattr(self, name, value)
            return

        if len(name) > 1:
            k = name[0]
            if hasattr(self, k):
                self.set_value(name[1:], value, getattr(self, k))
                return
            elif p:
                self.set_value(name[1:], value, p[k])
                return
            self.set_value(name[1:], value, self)
            return
        if p:
            p[name[0]] = value
            return
        setattr(self, name[0], value)

    def get(self, name, r=None):
        """return value of key name

        Args:
            name: key of parameter value

        Return:
            value of parameter identified by key
        """
        try:
            if isinstance(name, str):
                if r is not None:
                    return getattr(self, name, r)
                return getattr(self, name)
            if r and type(r) == dict:
                for k in name:
                    r = r.get(k)
                return r
            if len(name) > 1:
                if r:
                    return self.get(name[1:], getattr(r, name[0]))
                return self.get(name[1:], getattr(self, name[0]))
            return getattr(self, name[0])
        except KeyError as e:
            logger.error(e)
            raise MCerror(e)

    def __getitem__(self, name):
        return getattr(self, name)

    def __str__(self):
        "return string format of this object"
        return repr(self.__dict__)

    def __repr__(self):
        "representation of this object"
        return self.__str__()


class MachineModel(Model):
    """represents a machine model for a FE analysis

    Args:
      parameters: string or dict containing the model parameters. For example:
    ::

        {'lfe': 0.1,
        'poles': 4,
        'outer_diam': 0.13,
        'bore_diam': 0.07,
        'stator':{
        'num_slots': 12,
        'num_slots_gen': 3,
        ...
        },
        'magnet':{
            'material': 'M395',
            ..
        }
        }

        if parameters is string it is interpreted as the model name.
    """

    def __init__(self, parameters):
        super(self.__class__, self).__init__(parameters)
        name = 'DRAFT'
        if isinstance(parameters, str):
            name = parameters
        elif 'name' in parameters:
            name = parameters['name']
        # connect model even for complete model (see fsl connect_models)
        self.connect_full = True
        # must sanitize name to prevent femag complaints
        self.name = ''.join([n
                             for n in name.strip()
                             if n in MODELNAME_CHARS][:MAX_MODEL_NAME_LEN])
        try:
            self.external_rotor = (self.external_rotor == 1)
        except AttributeError:
            self.external_rotor = False
        self.move_inside = 1.0 if self.external_rotor else 0.0
        if 'exit_on_end' not in parameters:
            self.exit_on_end = 'false'
        if 'magnet' in parameters:
            if 'mcvkey_mshaft' in self.magnet:
                self.magnet['mcvkey_shaft'] = self.magnet['mcvkey_mshaft']
            for mcv in ('mcvkey_yoke', 'mcvkey_shaft'):
                if mcv not in self.magnet:
                    self.magnet[mcv] = 'dummy'
        if 'coord_system' in parameters:
            self.move_action = 1
        else:
            self.coord_system = 0
            self.move_action = 0

        try:
            self.set_num_slots_gen()
        except (AttributeError):
            pass
        try:
            if 'wdgtype' not in self.windings:
                self.windings['wdgtype'] = 'SYM'
        except AttributeError:
            pass

        try:
            self.commutator = self.windings['wdgtype'] == 'CMM'
        except AttributeError:
            self.commutator = False
        try:
            self.windings['cufilfact'] = self.windings['fillfac']
        except (KeyError, AttributeError):
            pass
        try:
            if 'leak_dist_wind' in self.windings:
                self.windings['leak_dist_wind']['wiredia'] = self.windings['dia_wire']
        except:
            pass

    def set_num_slots_gen(self):
        if 'num_slots_gen' not in self.stator:
            try:
                m = self.windings['num_phases']
            except (KeyError, AttributeError):
                m = 1

            if hasattr(self, 'magnet'):
                slotgen = [m*self.poles]
            else:
                slotgen = [self.poles]
            try:
                slotgen.append(self.stator['num_slots'])
            except KeyError:
                pass

            try:
                slotgen.append(self.rotor['num_slots'])
            except (AttributeError, KeyError):
                pass

            g = np.gcd.reduce(slotgen)
            if hasattr(self, 'magnet') and g > 1:
                g /= m
            if 'num_slots' in self.stator:
                Q1 = self.stator['num_slots']
                self.stator['num_slots_gen'] = Q1 // g
            else:
                try:
                    Q2 = self.rotor['num_slots']
                    self.rotor['num_slots_gen'] = Q2 // g
                except (AttributeError, KeyError):
                    pass

    def set_mcvkey_magnet(self, mcvkey):
        self.mcvkey_magnet = mcvkey

    def get_mcvkey_magnet(self):
        try:
            return self.mcvkey_magnet
        except AttributeError:
            return ''

    def _set_magnet_material(self, comp, magcurves, magnetmat):
        names = []
        missing = []
        magnet = 0
        mcv = 0
        try:
            if magnetmat:
                magnet = magnetmat.find(comp['material'])
            if magnet and 'mcvkey' in magnet:
                if magcurves:
                    mcv = magcurves.find(magnet['mcvkey'])
                if mcv:
                    logger.debug('magnet mcv %s', mcv)
                    comp['mcvkey_magnet'] = mcv
                    names.append((mcv, 1.0))
                else:
                    missing.append(magnet['mcvkey'])
                    logger.error('magnet mcv %s not found',
                                 magnet['mcvkey'])
        except KeyError:
            pass
        except AttributeError:
            if 'material' in comp:
                missing.append(comp['material'])
                logger.error('magnet material %s not found',
                             comp['material'])
        return names, missing

    def set_magcurves(self, magcurves, magnetmat={}):
        """set and return real names of magnetizing curve material

        Args:
            magcurves: :class: 'MagnetizingCurve' including
                              magnetizing curve materials

        Return:
            set of magnetizing curve names with fillfac that have to be created

        """
        names = []
        missing = []

        mcv = 0
        if 'stator' in self.__dict__:
            fillfac = self.stator.get('fillfac', 1.0)
            for k in ('mcvkey_yoke', 'mcvkey_teeth'):
                mcvname = k+'_name'
                if mcvname not in self.stator:
                    try:
                        if self.stator[k] != 'dummy':
                            if magcurves:
                                mcv = magcurves.find(self.stator[k])
                            if mcv:
                                logger.debug('stator mcv %s', mcv)
                                self.stator[k] = magcurves.fix_name(
                                    mcv, fillfac)
                                names.append((mcv, fillfac))
                                self.stator[mcvname] = mcv
                            else:
                                missing.append(self.stator[k])
                                logger.error('stator mcv %s not found',
                                             self.stator[k])
                    except KeyError:
                        pass

                else:  # mcvname in self.stator:
                    names.append((self.stator[mcvname], fillfac))
            n, m = self._set_magnet_material(self.stator, magcurves, magnetmat)
            names += n
            missing += m

        if 'magnet' in self.__dict__:
            rotor = self.magnet
            subregion = 'magnet'
        elif 'rotor' in self.__dict__:
            rotor = self.rotor
            subregion = 'rotor'
        else:
            # no magnet no rotor
            return set(names)

        fillfac = rotor.get('fillfac', 1.0)
        try:
            if (rotor['mcvkey_yoke'] != 'dummy' and
                    'mcvkey_yoke_name' not in rotor):
                if magcurves:
                    mcv = magcurves.find(rotor['mcvkey_yoke'])
                if mcv:
                    logger.debug('%s mcv %s',
                                 subregion, mcv)
                    rotor['mcvkey_yoke'] = magcurves.fix_name(
                        mcv, fillfac)
                    rotor['mcvkey_yoke_name'] = mcv
                    names.append((mcv, fillfac))
                else:
                    missing.append(rotor['mcvkey_yoke'])
                    logger.error('%s mcv %s not found',
                                 subregion,
                                 rotor['mcvkey_yoke'])
            elif 'mcvkey_yoke_name' in rotor:
                names.append((rotor['mcvkey_yoke_name'], fillfac))

        except KeyError:
            pass

        try:
            if (rotor['mcvkey_shaft'] != 'dummy' and
                    'mcvkey_shaft_name' not in rotor):
                if magcurves:
                    mcv = magcurves.find(rotor['mcvkey_shaft'])
                if mcv:
                    logger.debug('shaft mcv %s', mcv)
                    rotor['mcvkey_shaft'] = magcurves.fix_name(mcv)
                    self.stator['mcvkey_shaft_name'] = mcv
                    names.append((mcv, 1.0))
                else:
                    missing.append(rotor['mcvkey_shaft'])
                    logger.error(' shaft %s not found',
                                 rotor['mcvkey_shaft'])
            elif 'mcvkey_shaft_name' in rotor:
                names.append((self.stator['mcvkey_shaft_name'], 1.0))

        except KeyError:
            pass

        if 'magnet' in self.__dict__:
            n, m = self._set_magnet_material(self.magnet, magcurves, magnetmat)
            missing += m
            names += n
        if missing:
            raise MCerror("MC pars missing: {}".format(
                ', '.join(set(missing))))
        return set(names)

    def statortype(self):
        """return type of stator slot"""
        for k in self.stator:
            if isinstance(self.stator[k], dict):
                return k
        raise AttributeError("Missing stator slot model in {}".format(
            self.stator))

    def rotortype(self):
        """return type of rotor slot"""
        for k in self.rotor:
            if isinstance(self.rotor[k], dict):
                return k
        raise AttributeError("Missing rotor model in {}".format(self.magnet))

    def magnettype(self):
        """return type of magnet slot"""
        for k in self.magnet:
            if k != 'material' and isinstance(self.magnet[k], dict):
                return k
        raise AttributeError("Missing magnet model in {}".format(self.magnet))

    def is_complete(self):
        """check completeness of models"""
        if self.is_dxffile():
            return True
        try:
            self.statortype()
            if hasattr(self, 'rotor'):
                self.rotortype()
            else:
                self.magnettype()
            return True
        except AttributeError:
            return False

    def is_dxffile(self):
        if 'dxffile' in self.__dict__:
            if isinstance(self.dxffile, dict):
                return True
        return False

    def has_magnet(self):  # either rotor or magnet
        return hasattr(self, 'magnet')


class FeaModel(Model):
    def __init__(self, parameters):
        self.recsin = ''  # recalc mcv for dynamic simulation
        self.cufilfact = 0.45
        self.culength = 1.4
        self.wind_temp = 20
        self.slot_indul = 0.0
        self.skew_angle = 0.0
        self.num_skew_steps = 0
        self.num_par_wdgs = 1
        self.eval_force = 0
        self.optim_i_up = 0
        self.plots = []
        self.airgap_induc = []
        super(self.__class__, self).__init__(parameters)
        if parameters.get('calculationMode', '') == 'asyn_motor':
            self.recsin = 'flux'

    def get_num_cur_steps(self):
        """returns number of curSteps (used for progress calc)"""
        try:
            if self.calculationMode == 'psd_psq_fast':
                return self.round((self.maxiq-self.miniq)/self.delta_iq) + 1
            if self.calculationMode == 'ld_lq_fast':
                return self.num_cur_steps+1  # must include 0
        except:
            pass
        return 1
