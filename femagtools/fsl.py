"""
    femagtools.fsl
    ~~~~~~~~~~~~~~

    Creating FSL Scripts



"""
import logging
import mako
import mako.lookup
import os
import re
import sys
from femagtools.dxfsl.converter import convert
from . import __version__
logger = logging.getLogger(__name__)


class FslBuilderError(Exception):
    pass


class Builder:
    def __init__(self):
        self.lookup = mako.lookup.TemplateLookup(
            directories=[os.path.join(os.path.dirname(__file__), 'templates'),
                         os.path.join(os.getcwd(), '.')],
            disable_unicode=False, input_encoding='utf-8',
            output_encoding='utf-8',
            default_filters=['decode.utf8'])

        self.fsl_stator = False
        self.fsl_magnet = False

    def create_stator_model(self, model):
        mcv = ["mcvkey_yoke = '{}'"
               .format(model.stator.get('mcvkey_yoke', 'dummy')),
               "mcvkey_shaft = '{}'"
               .format(model.stator.get('mcvkey_shaft', 'dummy'))]
        if model.stator.get('mcvkey_teeth'):
            mcv.append("mcvkey_teeth = '{}'".format(
                model.stator['mcvkey_teeth']))

        return mcv + self.render_stator(model)

    def prepare_stator(self, model):
        templ = model.statortype()
        if templ != 'dxffile':
            return

        logger.info("Conv stator from %s",
                    model.stator['dxffile']['name'])
        params = {}
        params['split'] = model.stator[templ].get('split', False)
        params['show_plots'] = model.stator[templ].get('plot', False)
        params['write_fsl'] = True
        params['airgap'] = -1.0
        pos = 'in' if model.external_rotor else 'out'
        params['part'] = ('stator', pos)
        conv = convert(model.stator['dxffile']['name'], **params)

        model.stator['num_slots'] = conv.get('tot_num_slot')
        self.set_diameter_parameter(model, conv)
        if model.get('dy1'):
            model.set_value('outer_diam', model.get('dy1'))
            model.set_value('bore_diam', model.get('da1'))

        model.stator['dxf'] = dict(fsl=conv['fsl'])
        self.fsl_stator = True
        del model.stator[templ]

    def set_diameter_parameter(self, model, conv):
        for v in ('dy1', 'da1', 'da2', 'dy2'):
            try:
                model.set_value(v, conv[v])
            except KeyError:
                pass

    def prepare_diameter(self, model):
        dy1 = model.get('dy1', 0.0)
        if dy1:
            model.set_value('outer_diam', dy1 * 1e-3)
        da1 = model.get('da1', 0.0)
        if da1:
            model.set_value('bore_diam', da1 * 1e-3)
        dy2 = model.get('dy2', 0.0)
        if dy2:
            model.set_value('inner_diam', dy2 * 1e-3)
        da2 = model.get('da2', 0.0)
        if da1 and da2:
            model.set_value('airgap', abs(da1 - da2)/2/1e3)

    def render_stator(self, model):
        templ = model.statortype()
        if templ == 'dxf':
            return [
                'agndst = {}'.format(model.agndst*1e3),
                'ndt(agndst)'
            ] + model.stator['dxf']['fsl']
        if templ == 'statorFsl':
            #  obsolete
            if 'parameter' in model.stator['statorFsl']:
                return self.render_template(
                    model.stator['statorFsl']['content_template'],
                    model.stator['statorFsl']['parameter'])
            if isinstance(model.stator['statorFsl']
                          ['content_template'], str):
                with open(model.stator['statorFsl']
                          ['content_template']) as f:
                    templ = [l.strip() for l in f.readlines()]
            else:
                templ = model.stator['statorFsl']['content_template']
            return self.render_template(
                '\n'.join(templ),
                model.stator['statorFsl'])
        
        statmodel = model.stator.copy()
        statmodel.update(model.stator[templ])
        statmodel.update({
            'zeroangle': model.stator.get('zeroangle', 0),
            'rlength': model.stator.get('rlength', 1),
            'num_layers': model.windings.get('num_layers', 1)})
        
        fslcode = self.__render(statmodel, templ, stator=True)
        if fslcode:
            return (fslcode +
                    ['post_models("nodedistance", "ndst" )',
                     'agndst=ndst[1]*1e3'])

        logger.error('File {}.mako not found'.format(templ))
        return []

    def create_magnet_model(self, model):
        mcv = ["mcvkey_yoke = '{}'"
               .format(model.magnet.get('mcvkey_yoke', 'dummy')),
               "mcvkey_shaft = '{}'"
               .format(model.magnet.get('mcvkey_shaft', 'dummy'))]

        if 'magnetFsl' in model.magnet:
            #  obsolete
            if 'parameter' in model.magnet['magnetFsl']:
                return mcv + self.render_template(
                    model.magnet['magnetFsl']['content_template'],
                    model.magnet['magnetFsl']['parameter'])
            if isinstance(model.magnet['magnetFsl']
                          ['content_template'], str):
                with open(model.magnet['magnetFsl']
                          ['content_template']) as f:
                    templ = [l.strip() for l in f.readlines()]
            else:
                templ = model.magnet['magnetFsl']['content_template']

            return mcv + self.render_template(
                '\n'.join(templ),
                model.magnet['magnetFsl'])

        templ = model.magnettype()
        magmodel = model.magnet.copy()
        magmodel.update(model.magnet[templ])
        magmodel['mcvkey_magnet'] = model.get_mcvkey_magnet()
        if templ == 'dxf':
            return mcv + [
                u'xmag = {}',
                u'ymag = {}',
                u'mag_orient = {}',
                u'ndt(agndst)'] + model.magnet['dxf']['fsl']
            
        return mcv + self.render_rotor(magmodel, templ)

    def prepare_rotor(self, model):
        templ = model.magnettype()
        if templ != 'dxffile':
            return

        params = {}
        params['split'] = model.magnet[templ].get('split', False)
        params['show_plots'] = model.magnet[templ].get('plot', False)
        params['write_fsl'] = True
        params['airgap'] = -1.0
        pos = 'out' if model.external_rotor else 'in'
        params['part'] = ('rotor', pos)
        logger.info("Conv rotor from %s", templ + '.dxf')
        conv = convert(model.magnet[templ]['name'], **params)
        model.set_value('poles', int(conv.get('num_poles')))
        self.set_diameter_parameter(model, conv)
        if model.get('da2'):
            logger.info('da2 %f',  model.get('da2')/1e3)
            ag = (model.get('bore_diam') - model.get('da2')/1e3)/2
            model.set_value('airgap', ag)
            
        model.magnet['dxf'] = dict(fsl=conv['fsl'])
        self.fsl_magnet = True
        del model.magnet[templ]

    def render_rotor(self, magmodel, templ):
        fslcode = self.__render(magmodel, templ, magnet=True)
        if fslcode:
            return fslcode

        logger.error('File {}.fsl not found'.format(templ))
        return []

    def create_connect_models(self, model):
        """return connect_model if rotating machine"""
        if model.get('move_action') == 0:
            return ['pre_models("connect_models")']
        return []

    def create_open(self, model):
        return (['-- created by femagtools {}'.format(__version__), ''] +
                self.__render(model, 'open') +
                self.__render(model, 'basic_modpar'))

    def set_modpar(self, model):
        return self.__render(model, 'basic_modpar')

    def create_new_model(self, model):
        return (['-- created by femagtools {}'.format(__version__), ''] +
                self.__render(model, 'new_model'))

    def create_cu_losses(self, model):
        return self.__render(model.windings, 'cu_losses')

    def create_fe_losses(self, model):
        if model.get('ffactor', 0):
            return self.__render(model, 'FE-losses')
        return []
    
    def create_gen_winding(self, model):
        return self.__render(model, 'gen_winding')

    def create_model_with_dxf(self, model):
        dxfname = model.dxffile.get('name', None)
        if not dxfname:
            logger.error('Name of dxf-file expected')
            return []

        if dxfname.split('.')[-1] not in ('dxf'):
            dxfname += '.dxf'
        if not os.path.isfile(dxfname):
            logger.error('File {} not found'.format(dxfname))
            return []

        params = {}
        params['split'] = model.dxffile.get('split', False)
        params['show_plots'] = model.dxffile.get('plot', False)
        params['write_fsl'] = True
        params['airgap'] = model.dxffile.get('airgap', 0.0)

        conv = convert(dxfname, **params)

        model.set_value('poles', conv.get('num_poles'))
        model.set_value('outer_diam', conv.get('dy1') * 1e-3)
        model.set_value('bore_diam', conv.get('da1') * 1e-3)
        model.set_value('inner_diam', conv.get('dy2') * 1e-3)
        model.set_value('airgap', (conv.get('da1') - conv.get('da2'))/2/1e3)
        model.set_value('agndst', conv.get('agndst')*1e-3)

        if not hasattr(model, 'stator'):
            setattr(model, 'stator', {})
        model.stator['num_slots'] = conv.get('tot_num_slot')
        model.stator['num_slots_gen'] = conv.get('num_sl_gen')
        if 'fsl_stator' in conv:
            self.fsl_stator = True
            model.stator['dxf'] = dict(fsl=conv['fsl_stator'])
        if not hasattr(model, 'magnet'):
            setattr(model, 'magnet', {})
        if 'fsl_magnet' in conv:
            self.fsl_magnet = True
            model.magnet['dxf'] = dict(fsl=conv['fsl_magnet'])

        return self.create_new_model(model) + \
            self.create_cu_losses(model) + \
            self.create_stator_model(model) + \
            self.create_gen_winding(model) + \
            self.create_magnet(model) + \
            self.create_magnet_model(model) + \
            self.mesh_airgap(model) + \
            self.create_connect_models(model)

    def create_model(self, model):
        if model.is_dxffile():
            return self.create_model_with_dxf(model)

        if model.is_complete():
            self.prepare_stator(model)
            self.prepare_rotor(model)
            self.prepare_diameter(model)
            if self.fsl_stator:
                from femagtools.dxfsl.fslrenderer import agndst
                ag = model.get('airgap')
                model.set_value(
                    'agndst',
                    agndst(model.get('bore_diam'),
                           model.get('bore_diam') - 2*ag))

            model.set_num_slots_gen()

            return self.create_new_model(model) + \
                self.__render(model.windings, 'cu_losses') + \
                self.create_stator_model(model) + \
                self.__render(model, 'gen_winding') + \
                self.create_magnet(model) + \
                self.create_magnet_model(model) + \
                self.mesh_airgap(model) + \
                self.create_connect_models(model)

        return self.open_model(model)

    def open_model(self, model):
        return self.create_open(model)

    def load_model(self, model):
        return self.__render(model, 'open')

    def create_magnet(self, model, magnetMat=None):
        try:
            if magnetMat:
                if 'mcvkey' in magnetMat:
                    model.set_mcvkey_magnet(magnetMat['mcvkey'])
                return self.__render(magnetMat, 'magnet')
            return ['m.remanenc       =  1.2',
                    'm.relperm        =  1.05']
        except AttributeError:
            pass  # no magnet
        return []

    def create_common(self, model):
        return self.__render(model, 'common')

    def create_analysis(self, model, magnets=None, magnet_material=None):
        magnetMat = None
        magndata = []
        if magnets and magnet_material:
            magnetMat = magnets.find(magnet_material)
            if not magnetMat:
                raise FslBuilderError('magnet material {} not found'.format(
                    magnet_material))
            try:
                magnetMat['magntemp'] = model.magn_temp
            except AttributeError:
                pass
            magndata = self.create_magnet(model, magnetMat)

        airgap_induc = (self.create_airgap_induc()
                        if model.get('airgap_induc', 0) else [])
        
        if model.get('calculationMode') in ('cogg_calc',
                                            'ld_lq_fast',
                                            'pm_sym_loss',
                                            'torq_calc',
                                            'psd_psq_fast'):
            return (self.__render(model, model.get('calculationMode')) +
                    airgap_induc)
        
        return (self.__render(model, 'cu_losses') +
                magndata +
                self.__render(model, model.get('calculationMode')) +
                airgap_induc +
                self.__render(model, 'plots'))
#                ['save_model(cont)'])

    def create_airgap_induc(self):
            return self.__render(dict(), 'airgapinduc')

    def create_colorgrad(self, model):
            return self.__render(model, 'colorgrad')

    def mesh_airgap(self, model):
        if self.fsl_stator and self.fsl_magnet:
            return self.__render(model, 'mesh-airgap')
        else:
            return []

    def create(self, model, fea, magnets=None):
        "create model and analysis function"
        try:
            fea['lfe'] = model.get('lfe')
        except AttributeError:
            pass
        try:
            fea['move_action'] = model.get('move_action')
        except AttributeError:
            pass
        try:
            fea.update(model.windings)
        except AttributeError:
            pass
        
        if model.is_complete():
            logger.info("create new model and simulation")
            fslmodel = self.create_model(model)
            fea['pocfilename'] = (model.get('name') +
                                  '_' + str(model.get('poles')) +
                                  'p.poc')
            if not 'phi_start' in fea:
                fea['phi_start'] = 0.0
            if not 'range_phi' in fea:
                fea['range_phi'] = 720/model.get('poles')
            fe_losses = self.create_fe_losses(model)
            return fslmodel + fe_losses + self.create_analysis(
                fea, magnets, model.magnet.get('material', 0))
        
        logger.info("create open model and simulation")
        return self.open_model(model) + \
            self.create_analysis(fea)

    def __render(self, model, templ, stator=False, magnet=False):
        if templ.split('.')[-1] in ('fsl', 'mako'):
            try:
                template = self.lookup.get_template(templ)
                logger.info('use file {}'.format(templ))
                return template.render_unicode(model=model).split('\n')
            except mako.exceptions.TopLevelLookupException as ex:
                logger.error('File {} not found'.format(templ))
                sys.exit(1)

        try:
            template = self.lookup.get_template(templ+".mako")
            logger.debug('use template {}.mako'.format(templ))
        except mako.exceptions.TopLevelLookupException as ex:
            template = self.lookup.get_template(templ+".fsl")
            logger.debug('use FSL {}.fsl'.format(templ))
            if stator:
                self.fsl_stator = True
            if magnet:
                self.fsl_magnet = True

        return template.render_unicode(model=model).split('\n')

    def render_template(self, content_template, parameters):
        template = mako.template.Template(content_template)
        obj = {}
        if isinstance(parameters, list):
            for p in parameters:
                obj[p['key']] = p['value']
            return template.render_unicode(
                par=obj).split('\n')
        return template.render_unicode(
            model=parameters).split('\n')

    def read(self, fslfile):
        """extracts parameters from content and creates template"""
        parpat = re.compile(
            r'''(\w+)\s*= # key
            \s*(([+-]?\d+((?:\.\d+)?(?:[eE][+-]\d+)?))|
            (['"][^'"]*['"]))\s*-- # value
            \s*(.+)$''', re.X)
        content = []
        content_template = []
        parameter = []
        for line in fslfile:
            # add line to content
            content.append(line.strip())
            p = parpat.findall(line)
            if p:
                par = dict(key=p[0][0],
                           value=p[0][1].strip(),
                           comment=p[0][-1])
                parameter.append(par)
                content_template.append("{0} = {1} -- {2}".format(
                    par['key'],
                    '${par[\''+par['key']+'\']}',
                    par['comment']))
                continue

            # nothing do do => insert this line in template
            content_template.append(line.strip())

        return dict(
            content='\n'.join(content),
            type='fsl',
            content_template='\n'.join(content_template),
            parameter=parameter)
