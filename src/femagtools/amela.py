'''
    femagtools.amela
    ~~~~~~~~~~~~~~~~

    Calculate Magnet Losses with AMELA
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

def geometry_id(bndx, bndy):
    '''identify the magnet geometry''' 
    # caculate magnet area 
    xy = 0
    yx = 0
    for i in range(len(bndy)):
        if i < (len(bndy)-1):
            xy += bndx[i]*bndy[i+1]
            yx += bndy[i]*bndx[i+1]
        else:
            xy += bndx[i]*bndy[0]
            yx += bndy[i]*bndx[0]
    area = np.abs(yx-xy)*0.5

    x0, y0 = np.mean(bndx), np.mean(bndy)
    distances = [] 
    for i in range(len(bndx)): 
        for j in range(len(bndy)): 
            dist = np.sqrt((bndx[i] - bndx[j])**2 + (bndy[i] - bndy[j])**2)
            distances.append([dist, bndx[i], bndx[j], bndy[i], bndy[j]])
    
    distances.sort(reverse=True)
    xe = [distances[0][2], distances[2][1], distances[0][1], distances[2][2]]
    ye = [distances[0][4], distances[2][3], distances[0][3], distances[2][4]]

    # dimension 
    dim = np.zeros((3, 5))
    x1, y1 = xe[0], ye[0]

    for i in range(3): 
        x2, y2 = xe[i+1], ye[i+1]
        dim[i, 0] = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        dim[i, 1:3] = [x1, x2]
        dim[i, 3:5] = [y1, y2]

    dim = dim[np.lexsort((dim[:, 2], dim[:, 1], dim[:, 0]))]
    dx12, dy12 = dim[1, 1] - dim[1, 2],  dim[1, 3] - dim[1, 4]
    dx14, dy14 = dim[0, 1] - dim[0, 2],  dim[0, 3] - dim[0, 4]
    
    alp1 = np.arctan2(dy12, dx12)
    alp2 = np.arctan2(dy14, dx14)

    if alp1 < 0: alp1 += 2*np.pi
    if alp2 < 0: alp2 += 2*np.pi
    if alp2 < alp1: 
        alp2 += np.pi
        alpha = (alp1 + (alp2 - np.pi/2))/2
    else: 
        alpha = (alp1 + (alp2 - np.pi/2))/2 + np.pi
    
    wm = dim[1, 0]
    hm = area/wm 

    return dict(wm=wm, 
                hm=hm, 
                x0=x0, 
                y0=y0,
                area=area, 
                alpha=alpha)

def tf(b1, b2, alpha):
    '''Tranformation Matrix'''
    T = np.array([[cos(alpha), sin(alpha)], 
                [-sin(alpha), cos(alpha)]]) 
    if b1.ndim > 1: 
        r = T.dot(((b1.ravel()), (b2.ravel())))
        return [r[0, :].reshape(*b1.shape), 
                r[1, :].reshape(*b1.shape)]
    else: 
        return T.dot(((b1), (b2)))

def transform_coord(geometry, xcp, ycp): 
    '''transform from global coord to local coord'''
    # transformation 
    elcp = tf(b1=np.array(xcp)-geometry['x0'], 
              b2=np.array(ycp)-geometry['y0'], 
              alpha=geometry['alpha'])
    return dict(excpl=elcp[0, :]+geometry['wm']/2, 
                eycpl=elcp[1, :]+geometry['hm']/2, 
                excp=np.array(xcp), 
                eycp=np.array(ycp))

def transform_flux_denstiy(geometry, bx, by): 
    '''transform the magnet flux density to local coordinate system'''
    # transformation 
    bxy = tf(b1=bx, 
             b2=by, 
             alpha=geometry['alpha'])

    # remove DC component
    bxf = np.mean(bxy[0].T - np.mean(bxy[0],axis=1).T,axis=1)
    byf = np.mean(bxy[1].T - np.mean(bxy[1],axis=1).T,axis=1)               
     
    return dict(bxl=bxy[0],
                byl=bxy[1],
                bxf=bxf, 
                byf=byf 
                )

class Amela():
    '''Run Amela Calculation
    Parameters
    ----------
    workdir : str
        working directory of femag calculation
        (The directory that the nc file)
    amela_dir: str (optional)
        amela directory 
    magnet_data : dict
        calculation control
        -- name must be provided
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

    def __init__(self, workdir, magnet_data, amela_dir=None):

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

    def get_magnet_data(self, ibeta=None):
        '''Extract magnet data from nc file
        Parameters
        ----------
            None
        Returns
        ----------
        pm_data: list of magnet data
        '''
        nc_name = self.workdir / self.magn['name']
        r = femagtools.nc.read(nc_name)

        mag_spels = r.magnet_super_elements()
        spel_key = [i.key for i in mag_spels]
        spel_area = [float(i.area())*1e6 for i in mag_spels]

        pm_elem_key = []
        elem = []
        pm_node_key = [[] for i in range(len(spel_key))]
        bndkey = [[] for i in range(len(spel_key))]
        bx = [[] for i in range(len(spel_key))]
        by = [[] for i in range(len(spel_key))]
        xcp = [[] for i in range(len(spel_key))]
        ycp = [[] for i in range(len(spel_key))]
        bndx = [[] for i in range(len(spel_key))]
        bndy = [[] for i in range(len(spel_key))]
        # prepare data for ialh method
        wm = []
        hm = []
        alpha = []
        x0 = [] 
        y0 = []
        geometry = []
        elcp = []
        bl = []
        # conductivity and permeability of the magnets
        cond = 0
        mur = 0
        # read boundary nodes
        for k, i in enumerate(mag_spels):

            cond = i.conduc
            mur = 1/i.elements[0].reluc[0]
            logger.debug('Magnet: mur=%s, conductivity=%s', mur, cond)

            pm_elem_key.append([j.key for j in i.elements])
            elem.append(len(i.elements))
            for j in i.elements:
                for kk in j.vertices:
                    pm_node_key[k].append(kk.key - 1)

            for bnd in i.nodechains:
                for kk in bnd.nodes:
                    if len(bndkey[k]) > 0:
                        if kk.key != bndkey[k][-1]:
                            bndkey[k].append(kk.key)
                            bndx[k].append(kk.x*1e3)
                            bndy[k].append(kk.y*1e3)
                    else:
                        bndkey[k].append(kk.key)
                        bndx[k].append(kk.x*1e3)
                        bndy[k].append(kk.y*1e3)

                bndkey[k].pop(-1)
                bndx[k].pop(-1)
                bndy[k].pop(-1)
            geo = geometry_id(bndx=bndx[k], bndy=bndy[k])
            geometry.append(geo)
            # necessary?
            hm.append(geo['hm'])
            wm.append(geo['wm'])
            x0.append(geo['x0'])
            y0.append(geo['y0'])
            alpha.append(geo['alpha'])

        # default load angle (input beta I vs Up)
        if ibeta is not None:
            indx = ibeta
        else:
            num_cases = r.el_fe_induction_1.shape[3] - 1
            if 'loadcase' in self.magn:
                indx = self.magn['loadcase'] - 1
            else:
                indx = num_cases
            if indx == 3:
                indx = num_cases  # avoid error

        # stationary case, no rotation
        poles = 0
        try:
            poles = r.num_poles
        except:
            pass

        # read mesh and flux density
        for i in range(len(spel_key)):
            bx[i].append(np.zeros((elem[i], len(r.pos_el_fe_induction))))
            by[i].append(np.zeros((elem[i], len(r.pos_el_fe_induction))))
            for index, j in enumerate(mag_spels[i].elements):
                xcp[i].append(float(j.center[0]*1e3))
                ycp[i].append(float(j.center[1]*1e3))
                theta = np.arctan2(float(j.center[1]), float(j.center[0]))
                fd = r.flux_density(j, 0, indx)
                if poles == 0:
                    bx[i][0][index, :] = fd['bx']
                    by[i][0][index, :] = fd['by']
                else:
                    bx[i][0][index, :] = fd['bx']*np.cos(theta) - \
                        fd['by']*np.sin(theta)
                    by[i][0][index, :] = fd['bx']*np.sin(theta) + \
                        fd['by']*np.cos(theta)

            elcp.append(transform_coord(geometry[i], xcp[i], ycp[i]))
            bl.append(transform_flux_denstiy(geometry[i], bx[i][0], by[i][0]))


        if poles == 0:
            freq = self.magn.get('f', r.speed)
            time_vec = np.linspace(0, 1/freq, len(r.pos_el_fe_induction))
            pos = dict(time=time_vec.tolist(),
                       freq=freq,
                       t=float(1/freq))
            # reset num.poles
            poles = 1
        else:
            rpm = self.magn.get('speed', r.speed)
            ag_sim = r.pos_el_fe_induction[-1] - r.pos_el_fe_induction[0]
            pos = dict(phi=(r.pos_el_fe_induction*180/pi).tolist(),
                       speed=rpm,
                       t=float(60/rpm*ag_sim/360))  # TODO
        # prep dictionary for the loss calculation
        pm_data = []
        for i in range(len(spel_key)):
            pm_data.append(dict(name='pm_data_se' + str(spel_key[i]),
                                hm=self.magn.get('hm', hm[i]),
                                wm=self.magn.get('wm', wm[i]),
                                lm=self.magn.get('lm', r.arm_length*1e3),
                                alpha=alpha[i],
                                ls=r.arm_length*1e3,
                                sigma=float(self.magn.get('sigma', cond)),
                                mur=float(self.magn.get('mur', mur)),
                                loadcase=self.magn.get('loadcase', indx),
                                numpoles=poles,
                                nodes=dict(),
                                elements=dict(),
                                bndkeys=bndkey[i],
                                bndx=[float(c) for c in bndx[i]],
                                bndy=[float(c) for c in bndy[i]],
                                bl=bl[i],
                                elcp=elcp[i],
                                area=spel_area[i], 
                                spel_key=spel_key[i]))
            pm_data[i].update(pos)

        for k in range(len(pm_node_key)):
            for i, j in enumerate(pm_node_key[k]):
                pm_data[k]['nodes'][str(j + 1)] = dict(
                    key=int(j + 1),
                    x=float(r.nodes[j].x*1e3),
                    y=float(r.nodes[j].y*1e3),
                    Az=float(r.nodes[j].vpot[0])
                )
            for i, j in enumerate(pm_elem_key[k]):
                pm_data[k]['elements'][str(j)] = dict(
                    key=int(j + 1),
                    xcp=xcp[k][i],
                    ycp=ycp[k][i],
                    Bx=bx[k][0][i, :].tolist(),
                    By=by[k][0][i, :].tolist()
                )

        if len(mag_spels) / r.poles_sim > 1:
            return pm_data
        else:
            return [pm_data[0]]

    def get_magnet_data_all(self, num_op): 
        '''get all magnet data for all loadcases'''
        pm_data = []
        for i in num_op:
            pmd = self.get_magnet_data(ibeta=i)
            pm_data.append(pmd)
        return pm_data

    def export_json(self, pm_data):
        '''Export magnet data to json files
        Parameters
        ----------
        pm_data: list of magnet data
        Returns
        ----------
        None
        ----------
        '''
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

    def read_loss(self, pm_data):
        '''Read magnet losses
        Parameters
        ----------
            pm_data : dict
        Returns
        ----------
            losses : dict
        ----------
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

    def __call__(self, ialh=False):
        '''Run amela calculation
        Parameters
        ----------
            None
        Returns
        ----------
            losses : dict
        ----------
        '''
        # get magnet data
        r = self.get_magnet_data()
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
