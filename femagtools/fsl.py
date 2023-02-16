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
import math
import femagtools.windings
from femagtools.poc import Poc
from . import __version__
logger = logging.getLogger(__name__)


class FslBuilderError(Exception):
    pass


def cosys(model):
    if model.get('move_action', 0) == 0:
        return 'cosys("polar")'
    return 'cosys("cartes")'


class Builder:
    def __init__(self, templatedirs=[]):
        dirs = templatedirs
        if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
            # lookup up files in pyinstaller bundle
            logging.debug("Frozen!")
            dirs += [os.path.join(sys._MEIPASS, 'fsltemplates'),
                     os.path.join(os.getcwd(), '.')]
        else:
            dirs += [os.path.join(os.path.dirname(__file__), 'templates'),
                     os.path.join(os.getcwd(), '.')]
        self.lookup = mako.lookup.TemplateLookup(
            directories=dirs,
            input_encoding='utf-8',
            output_encoding='utf-8',
            default_filters=['decode.utf8'])

        self.fsl_stator = False
        self.fsl_magnet = False

    def create_wdg_def(self, model):
        name = 'winding'
        w = femagtools.windings.Winding(
            dict(
                Q=model.stator['num_slots'],
                p=model.poles//2,
                m=len(model.windings['wdgdef']),
                windings=model.windings['wdgdef']))
        self.copy_winding_file(name, w)
        return name

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
        if templ == 'mshfile':
            import femagtools.gmsh

            g = femagtools.gmsh.Gmsh(model.stator['mshfile']['name'])
            phi = g.get_section_angles()
            model.stator['num_slots'] = round(math.pi/(phi[1]-phi[0]))
            r = g.get_section_radius()
            model.set_value('outer_diam', 2*r[1])
            model.set_value('bore_diam', 2*r[0])
            model.stator['stator_msh'] = dict(
                name=model.stator['mshfile']['name']
            )
            for sr in g.get_subregions():
                model.stator[sr] = g.get_location(sr)
            for k in ('yoke', 'teeth', 'air'):
                model.stator[k] = model.stator['mshfile'].get(k, [])
            wdg = model.stator['mshfile'].get('wdg', '')
            if wdg:
                model.stator['wdg'] = model.stator[wdg]
            del model.stator['mshfile']

            self.fsl_stator = True
            return

        if templ == 'statorFsl':
            self.fsl_stator = True

        if templ != 'dxffile':
            return

        from femagtools.dxfsl.converter import convert
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
            elif model.stator['statorFsl'].get('content'):
                return (['agndst = {}'.format(model.get('agndst', 1e-3)*1e3),
                        'ndt(agndst)'] +
                        model.stator['statorFsl']['content'].split('\n'))
            if isinstance(model.stator['statorFsl']
                          ['content_template'], str):
                with open(model.stator['statorFsl']
                          ['content_template']) as f:
                    templ = [l.strip() for l in f.readlines()]
            else:
                templ = model.stator['statorFsl']['content_template']
            return self.render_template(
                '\n'.join(['agndst = {}'.format(model.agndst*1e3),
                           'ndt(agndst)'] + templ), model.stator['statorFsl'])

        statmodel = model.stator.copy()
        statmodel.update(model.stator[templ])
        k = 'zeroangle'
        if k not in statmodel:
            statmodel[k] = 0
        k = 'rlength'
        if k not in statmodel:
            statmodel[k] = 1
        k = 'num_layers'
        if k not in statmodel:
            statmodel[k] = 1
        fslcode = self.__render(statmodel, templ, stator=True)
        if fslcode:
            if self.fsl_stator:
                return (['agndst = {}'.format(model.agndst*1e3),
                         'alfa = 2*math.pi*m.num_sl_gen/m.tot_num_slot',
                         'num_agnodes = math.floor(m.fc_radius*alfa/agndst + 0.5)'] +
                        fslcode)
            if hasattr(model, 'num_agnodes'):
                return fslcode
            return (fslcode +
                    ['post_models("nodedistance", "ndst" )',
                     'agndst=ndst[1]*1e3'])

        logger.error('File {}.mako not found'.format(templ))
        return []

    def create_magnet_model(self, model):
        mcv = ["mcvkey_yoke  = '{}'"
               .format(model.magnet.get('mcvkey_yoke', 'dummy')),
               "mcvkey_shaft = '{}'"
               .format(model.magnet.get('mcvkey_shaft', 'dummy'))]

        if 'magnetFsl' in model.magnet:
            self.fsl_magnet = True
            #  obsolete
            if 'parameter' in model.magnet['magnetFsl']:
                return mcv + self.render_template(
                    model.magnet['magnetFsl']['content_template'],
                    model.magnet['magnetFsl']['parameter'])
            elif model.magnet['magnetFsl'].get('content'):
                return mcv + model.magnet['magnetFsl']['content'].split('\n')
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

    def create_rotor_model(self, model, condMat=[], ignore_material=False):
        mcv = ["mcvkey_yoke  = '{}'"
               .format(model.rotor.get('mcvkey_yoke', 'dummy')),
               "mcvkey_shaft = '{}'"
               .format(model.rotor.get('mcvkey_shaft', 'dummy'))]

        culosses = []
        templ = model.rotortype()
        rotmodel = model.rotor.copy()
        if not ignore_material:
            if 'conductivity' in rotmodel:
                culosses = self.create_cu_losses(
                    dict(cuconduct=rotmodel['conductivity'],
                         winding_inside=not model.external_rotor))
            elif 'material' in rotmodel:
                rotmodel['winding_inside'] = not model.external_rotor
                culosses = self.create_cu_losses(rotmodel, condMat)

        rotmodel.update(model.rotor[templ])
        rotmodel['is_rotor'] = True  # just in case for the template
        return mcv + culosses + self.render_rotor(rotmodel, templ)

    def create_rotor_winding(self, model):
        if hasattr(model, 'rotor') and model.rotortype() == 'rot_hsm':
            return self.render_rotor(model.rotor, 'rotor_winding')
        return []

    def prepare_magnet(self, model):
        templ = model.magnettype()
        if templ == 'mshfile':
            import femagtools.gmsh
            g = femagtools.gmsh.Gmsh(model.magnet['mshfile']['name'])
            r = g.get_section_radius()
            phi = g.get_section_angles()
            p = round(math.pi/(phi[1]-phi[0]))
            model.set_value('poles', p)
            model.set_value('inner_diam', 2*r[0])
            ag = (model.get('bore_diam') - 2*r[1])/2
            model.set_value('airgap', ag)
            model.magnet['rotor_msh'] = dict(
                name=model.magnet['mshfile']['name']
            )
            for sr in g.get_subregions():
                model.magnet[sr] = g.get_location(sr)

            for k in ('yoke', 'air'):
                model.magnet[k] = model.magnet['mshfile'].get(k, [])
            model.magnet['mag'] = dict(
                sreg=model.magnet['mshfile'].get('mag', []),
                axis=[g.get_axis_angle(s)
                      for s in model.magnet['mshfile'].get('mag', [])]
            )
            del model.magnet['mshfile']
            return

        if templ != 'dxffile':
            return

        from femagtools.dxfsl.converter import convert
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
        """return connect_model if rotating machine and incomplete model
        (Note: femag bug with connect model)"
        """
        if (model.get('move_action') == 0 and (
                model.connect_full or
                model.stator['num_slots'] > model.stator['num_slots_gen'])):
            return ['pre_models("connect_models")\n']
        return []

    def create_open(self, model):
        return (['-- created by femagtools {}'.format(__version__), ''] +
                self.__render(model, 'open') +
                ['global_unit("mm")', cosys(model)] +
                self.set_modpar(model) +
                self.create_fe_contr(model))

    def set_modpar(self, model):
        return self.__render(model, 'basic_modpar')

    def create_new_model(self, model):
        if model.get('num_agnodes', 0):
            tail = ['m.airgap = -2*ag/3']
        else:
            if isinstance(model.get('bore_diam', 0), list):
                tail = ['m.airgap   = 2*ag[2]/3']
            else:
                tail = ['m.airgap   = 2*ag/3']
        tail += [f"m.nodedist = {model.stator.get('nodedist',1)}"]

        return (['-- created by femagtools {}'.format(__version__), ''] +
                self.__render(model, 'new_model') +
                ['global_unit("mm")', cosys(model)] +
                self.set_modpar(model) +
                self.create_fe_contr(model) +
                tail)

    def create_fe_contr(self, model):
        return self.__render(model, 'fe-contr.mako')

    def create_cu_losses(self, windings, condMat=[], ignore_material=False):
        if 'material' in windings and not ignore_material:
            cond = 0
            if condMat:
                cond = condMat.find(windings['material'])
            if not cond:
                raise FslBuilderError(
                    'conductor material {} not found'.format(
                        windings['material']))
            windings['cuconduct'] = cond['elconduct']
        return self.__render(windings, 'cu_losses')

    def create_fe_losses(self, model):
        if any(model.get(k, 0) for k in ('ffactor', 'cw', 'ch', 'hyscoef',
                                         'edycof', 'indcof', 'fillfact',
                                         'basfreq', 'basind')):
            return self.__render(model, 'FE-losses')
        return []

    def create_gen_winding(self, model):
        genwdg = self.__render(model, 'gen_winding')
        k = list({'leak_dist_wind',
                  'leak_evol_wind',
                  'leak_tooth_wind'}.intersection(model.windings))
        if k:
            return (genwdg +
                    self.__render(model.windings[k[0]], k[0]) +
                    ['post_models("end_wind_leak","leak")',
                     'file_leak = io.open("end_wind_leak.dat","w")',
                     'file_leak:write(string.format("%g %g %g\\n", leak[1], leak[2], leak[3]))',
                     'file_leak:close()'])
        return genwdg

    def prepare_model_with_dxf(self, model):
        from femagtools.dxfsl.converter import convert
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
        params['nodedist'] = model.dxffile.get('nodedist', 1)

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

    def create_model(self, model, magnets=[], condMat=[], ignore_material=False):
        magnetMat = {}
        material = ''
        if not ignore_material:  # we are only looking for magnet material
            if hasattr(model, 'magnet'):
                material = model.magnet.get('material', 0)
            elif hasattr(model, 'stator') and not ignore_material:
                material = model.stator.get('material', 0)
        if magnets and material:
            magnetMat = magnets.find(material)
            if not magnetMat:
                raise FslBuilderError(
                    'magnet material {} not found'.format(
                        material))
            try:
                magnetMat['magntemp'] = model.magn_temp
            except AttributeError:
                magnetMat['magntemp'] = 20
        if model.is_complete():
            logger.info("create new model '%s'", model.name)
            if model.is_dxffile():
                self.prepare_model_with_dxf(model)
            else:
                self.prepare_stator(model)
                if hasattr(model, 'magnet'):
                    self.prepare_magnet(model)
                self.prepare_diameter(model)
                if self.fsl_stator:
                    from femagtools.dxfsl.fslrenderer import agndst
                    ag = model.get('airgap')
                    model.set_value(
                        'agndst',
                        agndst(model.get('bore_diam'),
                               model.get('bore_diam') - 2*ag,
                               model.stator.get('num_slots'),
                               model.get('poles'),
                               model.stator.get('nodedist') or 1.0))

                model.set_num_slots_gen()
            if hasattr(model, 'magnet'):
                if magnetMat:
                    model['magnet']['remanenc'] = magnetMat.get(
                        'remanenc', 1.2)
                    model['magnet']['relperm'] = magnetMat.get('relperm', 1.05)
                    model['magnet']['rlen'] = magnetMat.get('rlen', 1.0)
                rotor = (self.create_magnet(model) +
                         self.create_magnet_model(model))
                if magnetMat:
                    rotor += self.create_magnet(model, magnetMat)
            else:
                rotor = self.create_rotor_model(
                    model, condMat, ignore_material)
            try:
                windings = model.windings
            except:
                windings = {}
            windings['winding_inside'] = model.external_rotor
            if model.commutator:
                magdata = []
                if magnetMat:
                    magdata = self.create_magnet(model, magnetMat)

                return (self.create_new_model(model) +
                        self.create_cu_losses(windings, condMat, ignore_material) +
                        self.create_fe_losses(model) +
                        rotor +
                        ['m.remanenc    = 1.0',
                         'm.relperm     = 1.05'] +
                        self.create_stator_model(model) +
                        magdata +
                        self.create_gen_winding(model) +
                        self.mesh_airgap(model) +
                        self.create_connect_models(model) +
                        self.create_rotor_winding(model))
            if 'statorRing' in model.stator:
                return (self.create_new_model(model) +
                        self.create_cu_losses(windings, condMat, ignore_material) +
                        self.create_fe_losses(model) +
                        rotor +
                        self.create_stator_model(model) +
                        self.create_rotor_winding(model))

            return (self.create_new_model(model) +
                    self.create_cu_losses(windings, condMat, ignore_material) +
                    self.create_fe_losses(model) +
                    self.create_stator_model(model) +
                    self.create_gen_winding(model) +
                    rotor +
                    self.mesh_airgap(model) +
                    self.create_connect_models(model) +
                    self.create_rotor_winding(model))

        return (self.open_model(model) +
                self.create_fe_losses(model) +
                self.create_magnet(model, magnetMat))

    def open_model(self, model):
        return self.create_open(model)

    def load_model(self, model):
        return (self.__render(model, 'open') +
                ['global_unit("mm")', cosys(model)])

    def create_magnet(self, model, magnetMat={}):
        try:
            temp_prop = model.get('temp_prop',
                                  {k: magnetMat[k] for k in magnetMat})
            if temp_prop:
                logger.info("Setting magnet properties %s", magnetMat['name'])
                if hasattr(model, 'magnet'):
                    if 'rlen' in model.magnet:
                        temp_prop['rlen'] = model.magnet['rlen']
                elif hasattr(model, 'stator'):
                    if 'magn_rlen' in model.stator:
                        temp_prop['rlen'] = model.stator['magn_rlen']
                return self.__render(temp_prop, 'magnet-data')
            try:
                magnet = model.magnet
                rlen = magnet.get('rlen', 1)
            except AttributeError:
                magnet = model.stator  # commutator type?
                rlen = magnet.get('magn_rlen', 1)

            return ['m.remanenc       = {}'
                    .format(magnet.get('remanenc', 1.2)),
                    'm.relperm        = {}'
                    .format(magnet.get('relperm', 1.05)),
                    'm.rlen           = {}'
                    .format(100*rlen),
                    '']
        except:
            return []

    def create_analysis(self, sim):
        airgap_induc = (self.create_airgap_induc()
                        if sim.get('airgap_induc', 0) else [])
        felosses = self.create_fe_losses(sim)
        fslcalc = (self.__render(sim, sim.get('calculationMode')) +
                   airgap_induc)

        if sim.get('calculationMode') in ('cogg_calc',
                                          'ld_lq_fast',
                                          'pm_sym_loss',
                                          'torq_calc',
                                          'psd_psq_fast'):
            return felosses + fslcalc

        return (felosses + fslcalc +
                self.__render(sim, 'plots'))

    def create_shortcircuit(self, model):
        return self.__render(model, 'shortcircuit')

    def create_airgap_induc(self):
        return self.__render(dict(), 'airgapinduc')

    def create_colorgrad(self, model):
        return self.__render(model, 'colorgrad')

    def mesh_airgap(self, model):
        if ((self.fsl_stator and self.fsl_magnet) or
                model.get('num_agnodes', 0)):
            return self.__render(model, 'mesh-airgap')
        else:
            return []

    def create(self, model, sim, magnets=None, condMat=[]):
        "create model and analysis function"
        try:
            sim['lfe'] = model.get('lfe')
            num_poles = model.get('poles')
        except AttributeError:
            pass
        try:
            sim['move_action'] = model.get('move_action')
        except AttributeError:
            pass
        try:
            sim.update(model.windings)
            if 'num_poles' in model.windings:
                num_poles = model.windings['num_poles']
        except AttributeError:
            pass

        if 'poc' in sim:
            poc = sim['poc']
            poc.pole_pitch = 2*360/num_poles
            sim['pocfilename'] = poc.filename()
        elif 'pocfilename' not in sim:
            try:
                poc = Poc(2*360/num_poles)
                sim['poc'] = poc
                sim['pocfilename'] = poc.filename()
            except UnboundLocalError:
                logger.warning("unknown number of poles")
                pass

        if 'phi_start' not in sim:
            sim['phi_start'] = 0.0
        if 'range_phi' not in sim:
            try:
                sim['range_phi'] = 720/num_poles
            except UnboundLocalError:
                pass

        fslmodel = self.create_model(model, magnets, condMat)
        logger.info("create simulation '%s'", sim['calculationMode'])

        return (fslmodel + self.create_analysis(sim) +
                ['save_model("close")'])

    def __render(self, model, templ, stator=False, magnet=False):
        if templ.split('.')[-1] in ('fsl', 'mako'):
            try:
                template = self.lookup.get_template(templ)
                logger.debug('use file {}'.format(templ))
                return template.render_unicode(model=model).split('\n')
            except mako.exceptions.TopLevelLookupException:
                logger.error('File {} not found'.format(templ))
                sys.exit(1)

        try:
            template = self.lookup.get_template(templ+".mako")
            logger.debug('use template {}.mako'.format(templ))
        except mako.exceptions.TopLevelLookupException:
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
            r'''([\w\.]+)\s*= # key
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
