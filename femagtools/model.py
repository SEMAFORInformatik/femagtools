# -*- coding: utf-8 -*-
"""
    femagtools.model
    ~~~~~~~~~~~~~~~~

    Managing model parameters

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import logging

logger = logging.getLogger(__name__)


class MCerror(Exception):
    pass


class Model(object):
    def __init__(self, parameters):
        for k in parameters.keys():
            setattr(self, k, parameters[k])
            
    def set_value(self, name, value, p=None):
        """set value of parameter identified by name
        
        :param name: name of parameter
        :param value: value to be assigned to parameter
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

        :param name: key of parameter value
        :returns: value of parameter identified by key
        """        
        try:
            if isinstance(name, str):
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
     
    def __str__(self):
        "return string format of this object"
        return repr(self.__dict__)

    def __repr__(self):
        "representation of this object"
        return self.__str__()


class MachineModel(Model):
    """represents a machine model for a FE analysis
    
    :param parameters: dict containing the model parameters. For example:
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

    """
    def __init__(self, parameters):
        super(self.__class__, self).__init__(parameters)

        if 'name' in parameters:
            # must replace white space
            name = parameters['name'].strip().replace(' ', '_')
            for c in ('"', '(', ')'):
                name = name.replace(c, '')
            setattr(self, 'name', name)
        else:
            setattr(self, 'name', 'DRAFT')
        try:
            self.external_rotor = (self.external_rotor == 1)
        except:
            self.external_rotor = False
        self.move_inside = 1.0 if self.external_rotor else 0.0
        if 'magnet' in parameters:
            for mcv in ('mcvkey_yoke', 'mcvkey_mshaft'):
                if mcv not in self.magnet:
                    self.magnet[mcv] = 'dummy'
            
    def set_magcurves(self, magcurves):
        """set and return real names of magnetizing curve material

        :param magcurves: :class: 'MagnetizingCurve' including
                              magnetizing curve materials

        :returns: set of magnetizing curve names attached to this model

        """
        names = []
        missing = []
        if magcurves:
            if 'stator' in self.__dict__:
                try:
                    if self.stator['mcvkey_yoke'] != 'dummy':
                        mcv = magcurves.find(self.stator['mcvkey_yoke'])
                        if mcv:
                            logger.debug('stator mcv %s', mcv)
                            self.stator['mcvkey_yoke'] = mcv
                            names.append(mcv)
                        else:
                            missing.append(self.stator['mcvkey_yoke'])
                            logger.error('stator mcv %s not found',
                                         self.stator['mcvkey_yoke'])
                except KeyError:
                    pass
                
            if 'magnet' in self.__dict__:
                try:
                    if self.magnet['mcvkey_yoke'] != 'dummy':
                        mcv = magcurves.find(self.magnet['mcvkey_yoke'])
                        if mcv:
                            logger.debug('magnet mcv %s', mcv)
                            self.magnet['mcvkey_yoke'] = mcv
                            names.append(mcv)
                        else:
                            missing.append(self.magnet['mcvkey_yoke'])
                            logger.error('magnet mcv %s not found',
                                         self.magnet['mcvkey_yoke'])
                except KeyError:
                    pass

                try:
                    if self.magnet['mcvkey_mshaft'] != 'dummy':
                        mcv = magcurves.find(self.magnet['mcvkey_mshaft'])
                        if mcv:
                            logger.debug('mshaft mcv %s', mcv)
                            self.magnet['mcvkey_mshaft'] = mcv
                            names.append(mcv)
                        else:
                            missing.append(self.magnet['mcvkey_mshaft'])
                            logger.error('magnet shaft %s not found',
                                         self.magnet['mcvkey_mshaft'])
                except KeyError:
                    pass
        if missing:
            raise MCerror("MC pars missing: {}".format(
                ', '.join(set(missing))))
        return set(names)
    
    def statortype(self):
        """return type of stator slot"""
        for k in self.stator:
            if isinstance(self.stator[k], dict):
                return k
        raise MCerror("Missing stator slot model in {}".format(self.stator))
    
    def magnettype(self):
        """return type of magnet slot"""
        for k in self.magnet:
            if isinstance(self.magnet[k], dict):
                return k
        raise MCerror("Missing magnet model in {}".format(self.magnet))

    def is_complete(self):
        """check completeness of models"""
        try:
            self.statortype()
            self.magnettype()
            return True
        except:
            return False
        

class FeaModel(Model):
    def __init__(self, parameters):
        super(self.__class__, self).__init__(parameters)

    def __getitem__(self, k):
        return getattr(self, k)

if __name__ == '__main__':

    pars = dict(
         name = "PM 130 L4",
         lfe = 0.1,
         bore_diam = 0.07,
         outer_diam = 0.13,
         inner_diam = 0.04,
         airgap = 0.001,
         poles_gen = 1,
         poles = 4,
         stator = dict(
             zeroangle = 0.0,
             num_slots = 12,
             mcvkey_yoke = "3",
             num_slots_gen = 3,
             nodedist = 1.5,
             rlength = 1.0,
             statorRotor3 = dict(
                 slot_h1 = 0.002,
                 slot_h2 = 0.004,
                 middle_line = 0,
                 tooth_width = 0.009,
                 wedge_width2 = 0.0,
                 wedge_width1 = 0.0,
                 slot_top_sh = 0,
                 slot_r2 = 0.002,
                 slot_height = 0.02,
                 slot_r1 = 0.003,
                 slot_width = 0.003)),

         magnet = dict(
             mcvkey_yoke = "3",
             magn_len = 1.0,
             mcvkey_mshaft = "3",
             nodedist = 1.5,
             material = "1",
             
            magnetSector = dict(
                magn_num = 1,
                magn_width_pct = 0.8,
                magn_height = 0.004,
                magn_shape = 0.0,
                bridge_height = 0.0,
                magn_type = 1,
                condshaft_r = 0.02,
                magn_ori = 2,
                magn_rfe = 0.0,
                bridge_width = 0.0,
                magn_len = 1.0)),

          windings = dict(
               num_coils = 4.0,
               num_phases = 3,
               num_wires = 100,
               coil_span = 3.0,
               num_layers = 1) )

    model = Model( pars )
    print( model.get('lfe') )
    model.set_value(['name'],'Model2')
    model.set_value(['stator','yoke_diam'],0.133)
    print( model )
    model.set_value(['speed'],9000)
    print( model )
    print( model.get('name') )
#    print( model['name'] )
    print( model.get('u1') )

