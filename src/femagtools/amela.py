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

def T(alpha):
    return np.array([[np.cos(alpha), np.sin(alpha)],
                     [-np.sin(alpha), np.cos(alpha)]])

def transform_coord(geometry, xy):
    '''transform from global coord to local coord'''
    ecpl = T(geometry['alpha']).dot((xy-geometry['cxy']).T).T
    return dict(ecpl=(ecpl + (geometry['w']/2,
                              geometry['h']/2)).T,
                ecp=np.asarray(xy).T)

def tf(b1, b2, alpha):
    '''Transformation Matrix'''
    if b1.ndim > 1:
        r = T(alpha).dot(((b1.ravel()), (b2.ravel())))
        return [r[0, :].reshape(*b1.shape),
                r[1, :].reshape(*b1.shape)]
    else:
        return T(alpha).dot(((b1), (b2)))

def transform_flux_density(alpha, bxy):
    '''transform the magnet flux density to local coordinate system'''
    b = tf(b1=bxy[:, 0, :], b2=bxy[:, 1, :], alpha=alpha)

    # remove DC component
    bxf = np.mean(b[0].T - np.mean(b[0], axis=1).T, axis=1)
    byf = np.mean(b[1].T - np.mean(b[1], axis=1).T, axis=1)
    return {'bxyl': np.asarray(b),
            'bxyf': np.array([bxf, byf])}

def get_magnet_data(nc, ibeta=None) -> list:
    '''Extract magnet data from nc file

    Args:
        nc: nc object
        ibeta: load case

    Returns:
      pm_data: list of magnet data

    '''
    mag_spels = nc.magnet_super_elements()
    # prepare data for ialh method
    # conductivity and permeability of the magnets
    cond = 0
    mur = 0
    # read boundary nodes
    for se in mag_spels:
        cond = se.conduc
        if cond == 0:
            cond = 625000
            logging.info('Magnet conductivity equals 0, using 625000 S/m')
        mur = np.abs(1/se.elements[0].reluc[0])
        logging.debug('Magnet: mur=%s, conductivity=%s', mur, cond)

    # stationary case, no rotation
    poles = 0
    try:
        poles = nc.num_poles
    except AttributeError:
        pass

    if poles == 0:  # no rotation
        freq = nc.speed
        time_vec = np.linspace(0, 1/freq, len(nc.pos_el_fe_induction))
        pos = dict(time=time_vec.tolist(),
                   freq=freq,
                   t=float(1/freq))
        # reset num.poles
        poles = 1
    else:
        rpm = nc.speed
        ag_sim = nc.pos_el_fe_induction[-1] - nc.pos_el_fe_induction[0]
        pos = dict(phi=(nc.pos_el_fe_induction*180/np.pi).tolist(),
                   speed=rpm,
                   t=float(60/rpm*ag_sim/360))  # TODO
    # prep dictionary for the loss calculation
    pm_data = []
    for i, se in enumerate(mag_spels):
        ecp = [e.center for e in se.elements]
        geometry = se.get_rect_geom()

        bxy = []
        for e in se.elements:
            theta = np.arctan2(float(e.center[1]),
                               float(e.center[0]))
            fd = nc.flux_density(e, 0, ibeta)
            bxy.append(T(-theta).dot((fd['bx'], fd['by'])))
        #= np.moveaxis(bxy, 1, 0)
        pd = dict(name='pm_data_se' + str(se.key),
                  hm=geometry['h'],
                  wm=geometry['w'],
                  lm=nc.arm_length,
                  alpha=geometry['alpha'],
                  ls=nc.arm_length,
                  sigma=cond,
                  mur=mur,
                  loadcase=ibeta,
                  numpoles=poles,
                  nodes=dict(),
                  elements=dict(),
                  bl=transform_flux_density(geometry['alpha'],
                                            np.array(bxy)),
                  elcp=transform_coord(geometry, ecp),
                  area=se.area(),
                  spel_key=se.key)
        pd.update(pos)

        pm_data.append(pd)

    if len(mag_spels) / nc.poles_sim > 1:
        return pm_data
    else:
        return [pm_data[0]]

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

    def __init__(self, workdir: str, magnet_data: dict, amela_dir=None):

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
        return [get_magnet_data(r, ibeta=i) for i in num_op]

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

        r = get_magnet_data(nc, indx)
        # export to json
        self.export_json(r)
        # run amela
        calc_method = 'IALH' if ialh else '3DI'
        cmd = self.cmd + ['--calc', calc_method, self.magn['name'] + '/']
        log_file = self.amela_dir / 'amela.out'
        logger.info("Calculating magnet losses with AMELA ...")
        with log_file.open('w') as output:
            subprocess.run(cmd, stdout=output)
        return self.read_loss(r)
