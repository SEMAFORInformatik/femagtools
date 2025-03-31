'''Calculate Magnet Losses with AMELA

'''
import sys
import json
import pathlib
import logging
from pathlib import Path
import femagtools.nc
import numpy as np
from numpy import pi, sin, cos
import subprocess

#  set logging
logger = logging.getLogger(__name__)


class Amela():
    '''Run Amela Calculation

    Args:
      workdir : working directory of femag calculation
        (The directory that includes the nc file)
      amela_dir: str (optional) amela directory
      magnet_data : calculation control name must be provided
        dict(name='test',
        -- the following parameters are optional
             mur=1.05,
             sigma=625000,
             hm=3,
             wm=20,
             lm=30,
             speed=1000,
             nsegwid=0,
             nseglen=0)
    '''

    def __init__(self, workdir: str, magnet_data: dict,
                 amela_dir=None):

        self.magn = magnet_data
        self.workdir = pathlib.Path(workdir)
        if amela_dir is not None:
            self.amela_dir = pathlib.Path(amela_dir)
        else:
            self.amela_dir = self.workdir

        self.jsonfile = []
        if sys.platform == 'win32':
            self.cmd = [str(self.amela_dir / 'AMELA.BAT')]
        else:
            self.cmd = [str(self.amela_dir / 'AMELA')]
        # default batch
        self.cmd.append('--b')
        # append calc options
        if 'speed' in self.magn:
            self.cmd.append(f"--speed {self.magn['speed']}")
        if 'nsegwid' in self.magn:
            self.cmd.append(f"--nsegwid {self.magn['nsegwid']}")
        if 'nseglen' in self.magn:
            self.cmd.append(f"--nseglen {self.magn['nseglen']}")

    def get_magnet_data_all(self, num_op):
        '''get all magnet data for all loadcases'''
        nc_name = self.workdir / self.magn['name']
        r = femagtools.nc.read(nc_name)
        return [r.get_magnet_data(ibeta=i) for i in num_op]

    def export_json(self, pm_data: list):
        '''Export magnet data to json files

        Args:
          pm_data: list of magnet data

        '''
        return
        pm_dir = self.amela_dir / self.magn['name']
        pm_dir.mkdir(exist_ok=True)
        for i in pm_data:
            filename = (pm_dir / i['name']).with_suffix('.json')
            self.jsonfile.append(str(filename))  # for the future use
            # pop out non necessary data
            i.pop('bl')
            i.pop('elcp')
            with filename.open('w') as f:
                json.dump(i, f)
            logger.info('Exporting %s ...', i['name'])

    def read_loss(self, pm_data: dict) -> dict:
        '''Read magnet losses

        Args:
          pm_data : dict

        Returns:
          losses
        '''
        losses = {}
        for i in range(len(pm_data)):
            dirname = self.magn['name'] + '/' + pm_data[i]['name'] + \
                '/EClosses.tab'
            result_name = self.amela_dir / dirname
            with result_name.open() as f:
                data = np.loadtxt(f)
                total_loss = np.mean(data[0:-1, -1])
                losses[pm_data[i]['name']] = {"loss_data":data, "total_loss": total_loss}
                logger.info("Magnet losses in superelement %s is %.3f W",
                            pm_data[i]['name'].split('se')[-1], total_loss)
        return losses

    def __call__(self, ialh=False) -> dict:
        '''Run amela calculation

        Args:
          ialh: use method IALH if True else 3DI

        Returns:
            losses
        '''
        # get magnet data
        nc_name = self.workdir / self.magn['name']
        nc = femagtools.nc.read(nc_name)
        num_cases = nc.el_fe_induction_1.shape[3] - 1
        if 'loadcase' in self.magn:
            indx = self.magn['loadcase'] - 1
        else:
            indx = num_cases
        if indx == 3:
            indx = num_cases  # avoid error

        r = nc.get_magnet_data(indx)
        # export to json
        #self.export_json(r)
        # run amela
        calc_method = 'IALH' if ialh else '3DI'
        cmd = self.cmd + ['--calc', calc_method, self.magn['name'] + '/']
        log_file = self.amela_dir / 'amela.out'
        logger.info("Calculating magnet losses with AMELA ...")
        with log_file.open('w') as output:
            subprocess.run(cmd, stdout=output)
        return self.read_loss(r)
