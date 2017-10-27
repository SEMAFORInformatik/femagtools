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

logger = logging.getLogger(__name__)


class FslBuilderError(Exception):
    pass


class Builder:
    def __init__(self):
        self.lookup = mako.lookup.TemplateLookup(
            directories=[os.path.join(os.path.dirname(__file__), 'templates'),
                         os.path.join(os.getcwd(), 'templates')],
            disable_unicode=False, input_encoding='utf-8',
            output_encoding='utf-8',
            default_filters=['decode.utf8'])
    
    def create_stator_model(self, model):
        templ = model.statortype()
        return self.__render(model, templ)
    
    def create_magnet_model(self, model):
        mcv = ["mcvkey_yoke='{}'".format(
            model.magnet.get('mcvkey_yoke', 'dummy')),
               "mcvkey_shaft='{}'".format(
                   model.magnet.get('mcvkey_shaft', 'dummy'))]
        try:
            if 'magnetFsl' in model.magnet:
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
            return mcv + self.__render(model, templ)
        except AttributeError:
            pass  # no magnet
        return []

    def create_connect_models(self, model):
        """return connect_model if rotating machine"""
        if model.get('move_action') == 0:
            return ['pre_models("connect_models")']
        return []
    
    def create_open(self, model):
        return self.__render(model, 'open') + \
            self.__render(model, 'basic_modpar')

    def set_modpar(self, model):
        return self.__render(model, 'basic_modpar')

    def create_new_model(self, model):
        return self.__render(model, 'new_model')
    
    def create_model(self, model, magnets=None):
        magnetMat = None
        magndata = ['']
        if magnets and 'material' in model.magnet:
            magnetMat = magnets.find(model.magnet['material'])
            if not magnetMat:
                raise FslBuilderError('magnet material {} not found'.format(
                    model.magnet['material']))
            if not magnetMat.get('mcvkey', 0):
                magndata = ['pre_models("Magnet-data")', '']
        
        if model.is_complete():
            return self.create_new_model(model) + \
                self.__render(model.windings, 'cu_losses') + \
                self.create_stator_model(model) + \
                ['post_models("nodedistance", "ndst" )',
                 'agndst=ndst[1]*1e3'] + \
                self.__render(model, 'gen_winding') + \
                self.create_magnet(model, magnetMat) + \
                self.create_magnet_model(model) + \
                self.create_connect_models(model) + \
                magndata
        return self.open_model(model, magnets)

    def open_model(self, model, magnets=None):
        return self.create_open(model)
    
    def load_model(self, model, magnets=None):
        return self.__render(model, 'open')
    
    def create_magnet(self, model, magnetMat):
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
    
    def create_analysis(self, model):
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
                self.__render(model, model.get('calculationMode')) +
                airgap_induc +
                self.__render(model, 'plots') +
                ['save_model(cont)'])
            
    def create_airgap_induc(self):
            return self.__render(dict(), 'airgapinduc')
        
    def create_colorgrad(self, model):
            return self.__render(model, 'colorgrad')

    def create(self, model, fea, magnets=None):
        "create model and analysis function"
        try:
            fea['pocfilename'] = (model.get('name') +
                                  '_' + str(model.get('poles')) +
                                  'p.poc')
        except:
            pass
        try:
            fea['lfe'] = model.get('lfe')
        except:
            pass
        try:
            fea['move_action'] = model.get('move_action')
        except:
            pass
        try:
            fea.update(model.windings)
        except:
            pass
        if model.is_complete():
            return self.create_model(model, magnets) + \
                self.create_analysis(fea)
        return self.open_model(model, magnets) + \
            self.create_analysis(fea)
        
    def __render(self, model, templ):
        template = self.lookup.get_template(templ+".mako")
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

