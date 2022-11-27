#!/usr/bin/env python
#
import unittest
import femagtools.fsl
import femagtools.magnet
import copy
import re

modelpars = dict(
    name="PM 130 L4",
    outer_diam=0.13,
    bore_diam=0.07,
    lfe=0.1,
    poles=4,
    stator=dict(
        num_slots=12,
        mcvkey_yoke="3",
        num_slots_gen=3,
        nodedist=1.5,
        rlength=1.0),
    windings=dict(
        num_phases=3,
        num_layers=1,
        num_wires=4,
        coil_span=3))

feapars = dict(
    lfe=0.1,
    speed=50.0,
    current=10.0,
    nu_move_steps=49,
    num_cur_steps=5,
    angl_i_up=0,
    optim_i_up=0,
    wind_temp=60.0,
    magn_temp=60.0,
    eval_force=0,
    calc_fe_loss=1,
    cog_move_steps=90,
    num_layers=1,
    slot_indul=0,
    skew_angle=0.0,
    culength=1.4,
    num_par_wdgs=1,
    cufilfact=0.45,
    num_skew_steps=0)


class FslBuilderTest(unittest.TestCase):

    def setUp(self):
        self.m = copy.deepcopy(modelpars)
        self.builder = femagtools.fsl.Builder()

    def tearDown(self):
        self.m = None
        self.builder = None

    def test_stator1(self):
        self.m['stator']['stator1'] = dict(
            tooth_width=0.009,
            slot_rf1=0.002,
            tip_rh1=0.002,
            tip_rh2=0.002,
            slot_width=0.003)
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_stator_model(model)
        self.assertEqual(len(fsl), 27)

    def test_stator2(self):
        self.m['stator']['stator2'] = dict(
            slot_width=0.009,
            slot_t1=0.002,
            slot_t2=0.002,
            slot_t3=0.002,
            corner_width=0.002,
            slot_depth=0.003)
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_stator_model(model)
        self.assertEqual(len(fsl), 23)

    def test_stator3(self):
        self.m['stator']['statorRotor3'] = dict(
            slot_h1=0.002,
            slot_h2=0.004,
            middle_line=0,
            tooth_width=0.009,
            wedge_width2=0.0,
            wedge_width1=0.0,
            slot_top_sh=0,
            slot_r2=0.002,
            slot_height=0.02,
            slot_r1=0.003,
            slot_width=0.003)
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_stator_model(model)
        self.assertEqual(len(fsl), 40)

    def test_stator4(self):
        self.m['stator']['stator4'] = dict(
            slot_height=0.1,
            slot_h1=1e-3,
            slot_h2=0,
            slot_h3=2e-3,
            slot_h4=3e-4,
            slot_r1=11e-3,
            slot_width=22e-3,
            wedge_width1=111e-5,
            wedge_width2=222e-5,
            wedge_width3=333e-5)
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_stator_model(model)
        self.assertEqual(len(fsl), 39)

    def test_statorBG(self):
        self.m['stator']['statorBG'] = dict(
            yoke_diam_ins=0.0344,
            slottooth=0.0,
            tip_rad=0.0,
            middle_line=1,
            slot_h1=1e-3,
            slot_r1=0,
            slot_h3=2e-3,
            slot_r2=3e-4,
            tooth_width=3.2e-3,
            slot_width=22e-3)

        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_stator_model(model)
        self.assertEqual(len(fsl), 33)

    def test_magnetSector(self):
        self.m['magnet'] = dict(
            magnetSector=dict(
                magn_height=0.005,
                magn_width_pct=0.8,
                condshaft_r=0.0591,
                magn_rfe=0.0,
                magn_len=1.0,
                magn_shape=0.0,
                bridge_height=0.0,
                bridge_width=0.0,
                magn_ori=2,
                magn_type=1,
                magn_num=1))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        self.assertEqual(len(fsl), 30)

    def test_magnetIron(self):
        self.m['magnet'] = dict(
            magnetIron=dict(
                magn_height=0.005,
                magn_width=0.008,
                gap_ma_iron=0,
                air_triangle=5,
                iron_height=0.001,
                magn_rem=1.2,
                condshaft_r=0.0591,
                magn_ori=1,
                bridge_height=0,
                bridge_width=0,
                iron_shape=0))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        # print('\n'.join(fsl))
        self.assertEqual(len(fsl), 26)

    def test_magnetIron2(self):
        self.m['magnet'] = dict(
            magnetIron2=dict(
                magn_height=0.005,
                magn_width=0.008,
                gap_ma_iron=0,
                air_triangle=1,
                iron_height=0.001,
                magn_rem=1.2,
                condshaft_r=0.006,
                magn_ori=1,
                gap_ma_right=1e-3,
                gap_ma_left=2e-3,
                iron_shape=0))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        self.assertEqual(len(fsl), 28)

    def test_magnetIron3(self):
        self.m['magnet'] = dict(
            magnetIron3=dict(
                magn_height=0.005,
                magn_width=0.008,
                gap_ma_iron=0,
                iron_bfe=3e-3,
                magn_num=1,
                air_triangle=1,
                iron_height=0.001,
                condshaft_r=0.006,
                magn_ori=1,
                gap_ma_right=1e-3,
                gap_ma_left=2e-3,
                iron_shape=0))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        self.assertEqual(len(fsl), 27)

    def test_magnetIron4(self):
        self.m['magnet'] = dict(
            magnetIron4=dict(
                magn_height=0.005,
                magn_width=0.008,
                gap_ma_iron=0,
                iron_bfe=3e-3,
                magn_num=1,
                air_space_h=1e-3,
                iron_height=0.001,
                corner_r=0.006,
                magn_ori=1,
                magn_di_ra=1e-3,
                air_sp_ori=1,
                iron_shape=0))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        self.assertEqual(len(fsl), 27)

    def test_magnetIron5(self):
        self.m['magnet'] = dict(
            magnetIron5=dict(
                magn_height=0.005,
                magn_width=0.008,
                gap_ma_iron=0,
                iron_bfe=3e-3,
                magn_num=1,
                air_space_h=1e-3,
                iron_height=0.001,
                corner_r=0.006,
                air_space_b=1e-3,
                magn_di_ra=1e-3,
                air_sp_ori=1,
                iron_shape=0))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        self.assertEqual(len(fsl), 27)

    def test_magnetIronV(self):
        self.m['magnet'] = dict(
            magnetIronV=dict(
                magn_height=0.005,
                magn_width=0.008,
                gap_ma_iron=0,
                iron_hs=3e-3,
                magn_num=1,
                magn_rem=1.2,
                air_triangle=1,
                iron_height=0.001,
                condshaft_r=0.006,
                magn_angle=130,
                iron_shape=0))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        self.assertEqual(len(fsl), 27)

    def test_magnetFC2(self):
        self.m['magnet'] = dict(
            magnetFC2=dict(
                magn_height=0.005,
                magn_width=0.008,
                yoke_height=5e-3,
                iron_h1=3e-3,
                iron_h2=2e-3,
                iron_hp=2e-3,
                iron_b=2e-3,
                magn_num=1,
                iron_bfe=0.001,
                iron_bfo=0.001,
                iron_shape=0))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_magnet_model(model)
        self.assertEqual(len(fsl), 25)

    def test_rot_hsm(self):
        self.m['rotor'] = dict(
            rot_hsm=dict(
                gap_pol_shaft=1e-3,
                core_height=0.02,
                pole_height=0.016,
                pole_rad=0.042,
                core_width2=0.04,
                core_width1=0.04,
                pole_width_r=0.05,
                pole_width=0.052,
                slot_width=0.002,
                slot_height=0.002,
                damper_diam=0.004,
                damper_div=0.007
            ))
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_rotor_model(model, condMat=[])
        self.assertEqual(len(fsl), 33)

    def test_fe_losses(self):
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_fe_losses(model)
        self.assertFalse(fsl)
        self.m['ffactor'] = 1.1
        model = femagtools.MachineModel(self.m)
        fsl = self.builder.create_fe_losses(model)
        self.assertEqual(len(fsl), 21)
        ffact = [float(f.split('=')[1])
                 for f in fsl if f.startswith('m.ffactor')][0]
        self.assertEqual(ffact, self.m['ffactor'])
        feloss = [f.split('"')[1]
                  for f in fsl if f.find('pre_models') > 0][0]
        self.assertEqual(feloss, 'FE-Losses-1')

    def test_run_models(self):
        feapars['calculationMode'] = "cogg_calc"
        fsl = self.builder.create_analysis(feapars)
        self.assertEqual(len(fsl), 26)

        feapars['calculationMode'] = "pm_sym_fast"
        fsl = self.builder.create_analysis(feapars)
        self.assertEqual(len(fsl), 33)

        feapars['calculationMode'] = "mult_cal_fast"
        fsl = self.builder.create_analysis(feapars)
        self.assertEqual(len(fsl), 31)

        feapars['calculationMode'] = "torq_calc"
        fsl = self.builder.create_analysis(feapars)
        self.assertEqual(len(fsl), 26)

    def test_run_existing_model(self):
        model = femagtools.MachineModel('data/magnsec')
        feapars['calculationMode'] = "cogg_calc"
        fsl = self.builder.create(model, feapars)
        self.assertEqual(len(fsl), 65)

    def test_create_plots(self):
        pars = copy.deepcopy(feapars)
        pars['calculationMode'] = "pm_sym_fast"
        pars['plots'] = ['field_lines', 'Babs']
        fsl = self.builder.create_analysis(pars)
        field_lines = re.findall(r'field_lines\(([^)]*)\)', ''.join(fsl))
        self.assertEqual(len(field_lines), 1)
        self.assertEqual(int(field_lines[0].split(',')[-1]), 20)

        colorgrad = re.findall(r'color_gradation\(([^)]*)\)', ''.join(fsl))
        self.assertEqual(len(field_lines), 1)
        min, max = [int(l) for l in colorgrad[0].split(',')[4:6]]
        self.assertEqual(min, 0)
        self.assertEqual(max, 0)

        pars['plots'] = [('field_lines', 10), ('Babs', 0.0, 2.0)]
        fsl = self.builder.create_analysis(pars)
        field_lines = re.findall(r'field_lines\(([^)]*)\)', ''.join(fsl))
        self.assertEqual(len(field_lines), 1)
        self.assertEqual(int(field_lines[0].split(',')[-1]), 10)

        colorgrad = re.findall(r'color_gradation\(([^)]*)\)', ''.join(fsl))
        self.assertEqual(len(field_lines), 1)
        min, max = [float(l) for l in colorgrad[0].split(',')[4:6]]
        self.assertEqual(min, 0.0)
        self.assertEqual(max, 2.0)

    def test_readfsl(self):
        content = [
            'dshaft = 360 --shaft diameter',
            'hm  = 38 -- magnet height',
            'bm = 122 -- magnet width',
            'ws = 10  -- slot width',
            'lfe = 224',
            '-- calculate slot height, angle and pole pairs',
            'hs = (da2-dy2)/2 - bm   ',
            'alpha = math.pi/p/2     -- slot angle',
            'p   = m.num_poles/2',
            'x = {}',
            'y = {}',
            '-- Berechnung der Koordinaten',
            'x[1],y[1] = pr2c(dy2/2, 0)',
            'x[2],y[2] = pr2c(da2/2, 0)',
            'x[3],y[3] = pr2c(da2/2, alpha - math.atan2(ws,(da2/2)))',
            'x[4],y[4] = pr2c(da2/2-hs, alpha - math.atan2(ws,(da2/2 - hs)))',
            'nc_line(x[1], y[1], x[2], y[2], 0)']
        result = self.builder.read(content)
        self.assertEqual(len(result['parameter']), 4)
        for p in result['parameter']:
            self.assertTrue(p['key'] in ['dshaft', 'hm', 'bm', 'ws'])

    def test_gen_winding(self):
        model = femagtools.MachineModel(self.m)

        fsl = self.builder.create_gen_winding(model)
        self.assertEqual(len(fsl), 20)

        model.windings['leak_dist_wind'] = dict(
            perimrad=1,
            vbendrad=1,
            endheight=1,
            wiredia=1)
        fsl = self.builder.create_gen_winding(model)
        self.assertEqual(len(fsl), 33)

        model.windings.pop('leak_dist_wind')
        model.windings['leak_evol_wind'] = dict(
            evol1rad=1,
            evol2rad=1,
            botlevel=1,
            toplevel=1,
            endheight=1,
            evolbend=1,
            wiredia=1)
        fsl = self.builder.create_gen_winding(model)
        self.assertEqual(len(fsl), 38)

        model.windings.pop('leak_evol_wind')
        model.windings['leak_tooth_wind'] = dict(
            endheight=1,
            bendrad=1,
            wiredia=1)
        fsl = self.builder.create_gen_winding(model)
        self.assertEqual(len(fsl), 34)

    def test_create_model_with_magnet_material(self):
        magnetmat = [dict(
            name='M45',
            rlen=0.9,
            remanenc=1.1,
            relperm=1.04,
            spmaweight=7.4,
            temcoefbr=-0.0015,
            temcoefhc=-0.0013,
            magncond=625000.0
        )]

        machine = dict(
            name="PM 886 32",
            lfe=0.224,
            poles=32,
            outer_diam=0.886,
            bore_diam=0.76,
            inner_diam=0.4956,
            airgap=0.007,
            external_rotor=1,

            stator=dict(
                num_slots=120,
                rlength=1.0,
                stator4=dict(
                    slot_height=0.035,
                    slot_h1=0.002,
                    slot_h2=0.0,
                    slot_h3=0.004,
                    slot_h4=0.0,
                    slot_width=0.01,
                    slot_r1=0.0,
                    wedge_width1=0.01,
                    wedge_width2=0.01,
                    wedge_width3=0.01)
            ),

            magnet=dict(
                material='M45',
                magnetSector=dict(
                    magn_num=1,
                    magn_height=0.014,
                    magn_width_pct=0.85,
                    condshaft_r=0.0,
                    magn_rfe=0.0,
                    magn_len=1.0,
                    magn_shape=0.0,
                    bridge_height=0.0,
                    bridge_width=0.0,
                    magn_ori=1,
                    magn_type=1
                )
            ),

            windings=dict(
                material='Cu',
                num_phases=3,
                num_wires=5,
                coil_span=1,
                num_layers=2)
        )
        model = femagtools.MachineModel(machine)
        magnets = femagtools.magnet.Magnet(magnetmat)
        condMat = femagtools.magnet.Magnet([dict(name='Cu', elconduct=56e6)])
        fsl = self.builder.create_model(model, magnets, condMat)
        self.assertEqual(len(fsl), 185)
        brem = [l.strip() for l in fsl
                if l.split('=')[0].strip() == 'm.remanenc']
        self.assertEqual(brem[-1].split('=')[-1].strip(),
                         str(magnetmat[0]['remanenc']))
        rlen = [l.strip() for l in fsl
                if l.split('=')[0].strip() == 'm.rlen']
        self.assertEqual(rlen[0].split('=')[-1].strip(),
                         str(100*magnetmat[0]['rlen']))


if __name__ == '__main__':
    unittest.main()
