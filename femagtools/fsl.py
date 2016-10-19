"""
    femagtools.femag
    ~~~~~~~~~~~~~~~~

    Creating FSL Scripts

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
import logging
import mako
import mako.lookup
import os
import re
import math

logger = logging.getLogger(__name__)


class FslBuilderError(Exception):
    pass


class Builder:
    def __init__(self):
        self.lookup = mako.lookup.TemplateLookup(
            directories=[os.path.join(os.path.dirname(__file__), 'templates')],
            disable_unicode=False, input_encoding='utf-8',
            output_encoding='utf-8',
            default_filters=['decode.utf8'])
    
    def create_stator_model(self, model):
        templ = model.statortype()
        return self.__render(model, templ)
    
    def create_magnet_model(self, model):
        if 'magnetFsl' in model.magnet:
            if 'parameter' in model.magnet['magnetFsl']:
                return self.render_template(
                    model.magnet['magnetFsl']['content_template'],
                    model.magnet['magnetFsl']['parameter'])
            if isinstance(model.magnet['magnetFsl']['content_template'], str):
                with open(model.magnet['magnetFsl']['content_template']) as f:
                    templ = [l.strip() for l in f.readlines()]
            else:
                templ = model.magnet['magnetFsl']['content_template']
            logger.info(model.magnet['magnetFsl'])
            return self.render_template(
                '\n'.join(templ),
                model.magnet['magnetFsl'])
            
        templ = model.magnettype()
        return self.__render(model, templ)

    def create_open(self, model):
        return self.__render(model, 'open')
    
    def create_new_model(self, model):
        return self.__render(model, 'new_model')
    
    def create_model(self, model, magnets=None):
        return self.create_new_model(model) + \
            self.__render(model.windings, 'cu_losses') + \
            self.create_stator_model(model) + \
            ['post_models("nodedistance", "ndst" )',
             'agndst=ndst[1]*1e3'] + \
            self.__render(model, 'gen_winding') + \
            self.create_magnet(model, magnets) + \
            ["mcvkey_yoke='" + model.magnet['mcvkey_yoke'] + "'",
             "mcvkey_shaft='" + model.magnet['mcvkey_mshaft'] + "'"] + \
            self.create_magnet_model(model) + \
            ['pre_models("connect_models")\n']
    
    def open_model(self, model, magnets=None):
        return self.create_open(model) + \
            self.create_magnet(model, magnets) + \
            self.__render(model.windings, 'cu_losses')
    
    def create_magnet(self, model, magnets):
        if magnets and 'material' in model.magnet:
            magnet = magnets.find(model.magnet['material'])
            if magnet:
                return self.__render(magnet, 'magnet')
            raise FslBuilderError('magnet material {} not found'.format(
                model.magnet['material']))
        return [' m.remanenc       =  1.2',
                ' m.relperm        =  1.05']
    
    def create_common(self, model):
        return self.__render(model, 'common')
    
    def create_analysis(self, model):
        if model['calculationMode'] in ('cogg_calc',
                                        'ld_lq_fast',
                                        'pm_sym_loss',
                                        'psd_psq_fast'):
            return self.__render(model, model['calculationMode'])
        return self.__render(model, 'cu_losses') + \
            self.__render(model, model['calculationMode'])
            
    def create_airgap_induc(self, model):
            return self.__render(model, 'airgapinduc')
        
    def create_colorgrad(self, model):
            return self.__render(model, 'colorgrad')

    def create(self, model, fea, magnets=None):
        "create model and analysis function"
        fea['lfe'] = model.get('lfe')
        fea.update(model.windings)
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
    
    def createNtib(self, speed, current, beta, r1):
        """return Ntib info"""
        return ['Speed[1/min], Current[A], beta[Degr], R_STATOR[Ohm], n-p',
                'Number of values:      {0}'.format(len(speed))] + \
                ['{0:12.1f}{1:12.3f}{2:12.3f}{3:12f}{4:12f}'.format(
                    60*n, math.sqrt(2.0)*i1,
                    b,
                    self.machine['r1'], 3.0)
                 for n, i1, b in zip(speed, current, beta)]

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

