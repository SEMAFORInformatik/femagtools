'''Calculate Magnet Losses with IALH Method

'''
__author__ = 'Max Hullmann, Dapu Zhang, Ivan Solc'

import logging
import warnings
import numpy as np
from numpy import sinh, sin, cosh, cos, pi
from scipy.interpolate import RBFInterpolator

MUR0 = 4*np.pi*10**-7
#  set logging
logger = logging.getLogger(__name__)

class FrequencyDomain:
    def __init__(self, s):
        self.cmplx_amp = s['cmplx_amp']
        self.amp = s['amp']
        self.freq = s['freq']
        self.phase = s['phase']
        self.order = s['order']


def fd(s):
    f = FrequencyDomain(s)
    return f

def fft(s, fs):
    '''Calculate 1D FFT
    Parameters
    ----------
    s: time signal
    fs: sampling frequency
    Returns
    ----------
    cmplx_amp: complex amplitude
    amp:       amplitude (abs)
    freq:      frequency
    phase:     phase offset
    order:     order
    '''
    l = len(s)
    ll = l//2+1
    mag = np.fft.fft(s)

    tmp = 2*np.abs(mag)/l
    amp = tmp[0:ll]
    amp[0] = amp[0]/2
    amp[-1] = amp[-1]/2

    idx_max = np.argmax(amp[1:])+1
    freq = np.array([i*fs/l for i in range(ll)])

    return fd(dict(cmplx_amp=mag,
                   amp=amp,
                   freq=freq,
                   phase=np.angle(mag[0:ll]),
                   order=freq/freq[idx_max]))

def dfac(x):
    return 6/(x)**3*(sinh(x)-sin(x))/ \
                (cosh(x)+cos(x))

def nptns(x, x1, nx):
    ne = len(x)
    dx = 0.5/nx
    exy = np.sort(x1)
    ny = 2
    for i in range(ne-1):
        if 2*(exy[i+1] - exy[i])/(exy[i+1] + exy[i]) < dx:
            ny = ny + 1
        else:
            break
    return ny

def ngrid(wm, hm, elxy):
    '''calculate number of grid points for the interpolation'''
    be = np.sqrt(wm*hm/elxy['ecp'].shape[0])
    nx = np.around(wm/be) + 1
    ny = np.around(hm/be) + 1

    ny_new = nptns(elxy['excp'], elxy['excpl'], nx)
    nx_new = nptns(elxy['eycp'], elxy['eycpl'], ny)

    if np.abs(nx-nx_new)/nx < 0.5 and np.abs(ny-ny_new)/ny < 0.5:
        nx = nx_new
        ny = ny_new

    return [int(nx), int(ny)]

def binterp(x, y, xq, yq, b):
    '''interpolate flux density with Rbf interpolator'''
    f = RBFInterpolator(np.vstack((x, y)).T, b)
    inp = f(np.vstack((xq, yq)).T)
    return inp.reshape(len(np.unique(yq)), len(np.unique(xq)), -1)

def binterp_ialh2(x, y, xq, yq, b):
    '''interpolate flux density with Rbf interpolator'''
    yy = np.vstack((x, y)).T
    logger.debug("shape y %s b %s", yy.shape, np.shape(b))
    f = RBFInterpolator(yy, b, kernel='thin_plate_spline')
    inp = f(np.array([[i, j] for i, j in zip(xq, yq)]))
    return inp.reshape(len(np.unique(xq)), -1)

def Segmentation(wm, hm, lm, elxy, nsegx, nsegy, nsegz):
    ''' Creates new magnet elements' center points' x,y coordinates based on the number of segments and the original data resolution

    Inputs:  total magnet dimensions (width - wm, height - hm, axial length - lm)
             elxy - dictionary with global and local reference frame x,y coordinates of elements centerpoints
             number of magnet segments : nsegx - number of segments in x (wm), nsegy - segments in y = 1 (hm), nsegz - segments in axial direction (lm)
    Returns: number of elements in x direction per magnet segment - nx_new
             number of elements in y direction per magnet segment - ny new
             number of elements in z direction per magnet segment - nz new
             x,y coordinates of each element's centerpoint in local reference frame, for entire magnet - xx, yy
    '''
    # Default nx,ny,nz without considering the segmentation
    be = np.sqrt(wm*hm/np.shape(elxy['ecp'])[1])  #square elements
    nx_new = int(np.around(wm/be))
    ny_new = int(np.around(hm/be))
    nz_new = int(np.around(lm/be))

    if (nsegx > 1) or (nsegy > 1):
        # limits the number of elements per segment to 10 to avoid unnecessarily detailed resolution
        nx_new = int(max(10,int(np.around(nx_new*np.sqrt(1/nsegx)))))
        ny_new = int(max(10,int(np.around(ny_new*np.sqrt(1/nsegy)))))
        nz_new = int(max(10,int(np.around(nz_new*np.sqrt(1/nsegz)))))

    wms = wm/nsegx
    hms = hm/nsegy
    x0 = 0 # offset for ecpl
    y0 = 0 # offset for ecpl

    segxcpl = np.linspace(wms/2/nx_new, wms - wms/nx_new/2, nx_new) + x0 # x center points of wanted elements distribution, local ref frame, 1st segment
    segycpl = np.linspace(hms/2/ny_new, hms - hms/ny_new/2, ny_new) + y0 # y center points of wanted elements distribution, local ref frame, 1st segment

    x,y = np.meshgrid(segycpl, segxcpl)
    xx = np.zeros((x.shape[0]*nsegx, x.shape[1]*nsegy))  # new centerpoint coordinates for entire magnet
    yy = np.zeros_like(xx)
    a=np.zeros((x.shape[0]*nsegx))
    b=np.zeros((y.shape[1]*nsegy))

    for ii in range(nsegx):
        a[ii*nx_new: (ii+1)*nx_new] = segxcpl + wms*ii # generation of x coordinates of elements for each magnet segment
    for jj in range(ny_new*nsegy):
        xx[:,jj] = np.atleast_2d(a)
    for jj in range(nsegy):
        b[jj*ny_new: (jj+1)*ny_new] = segycpl + hms*jj # generation of y coordinates of elements for each magnet segment
    for ii in range(nx_new*nsegx):
        yy[ii,:] = np.atleast_2d(b)

    return nx_new, ny_new, nz_new, xx, yy

class MagnLoss:
    '''Calculate Magnet Losses with IALH Methode
    Parameters
    ----------
    nc: (object) nc/isa7
    ibeta: load cases [0, 1, 2]
    magnet_data: array of pm data (see nc/isa7 get_magnet_data)
    '''
    def __init__(self, **kwargs):
        #super().__init__(workdir, magnet_data=dict(name=modelname))
        if 'nc' in kwargs:
            ibeta = kwargs.get('ibeta', [0])
            nc = kwargs['nc']
            self.pm = [nc.get_magnet_data(ibeta=i) for i in ibeta]
        elif 'magnet_data' in kwargs:
            self.pm = kwargs['magnet_data']
        self.speed = kwargs.get('speed', self.pm[-1][-1]['speed'])
        logger.info("Speed %f", self.speed)
        try:  # move action rotation
            self.theta = self.pm[-1][-1]['phi'] # rotor pos
            self.lt = len(self.theta)
            self.tgrid = 60/self.speed*(self.theta[-1]-self.theta[0])/360
        except: # move action linear
            self.displ = self.pm[-1][-1]['displ'] # rotor pos
            self.lt = len(self.displ)
            self.tgrid = (self.displ[-1]-self.displ[0])/self.speed

        self.ls = self.pm[-1][-1]['ls']
        self.th_loss = []
        try:
            self.numpoles = self.pm[-1][-1]['numpoles']
        except KeyError:
            self.numpoles = 1

        try:
            self.mur = kwargs.get('mur', self.pm[-1][-1]['mur'])
        except KeyError:
            self.mur = 1.05

        self.sigma = kwargs.get('sigma', self.pm[-1][-1]['sigma'])
        self.symmetry = kwargs.get('symmetry', True)
        self.is_x = False
        self.segx = kwargs.get('segx', [1])

        try:
            self.lm = self.pm[-1][-1]['lm']
        except AttributeError:
            self.lm = 0

        self.segz = kwargs.get('segz', [0])
        # Determine the number of segments in z direction
        for i in range(len(self.segz)):
            if self.segz[i] > 0:
                self.lm = self.ls/self.segz[i]
            elif self.lm > 0:
                self.segz[i] = np.around(self.ls/self.lm)
            else:
                self.lm = self.ls
                self.segz[i] = 1

    def skin_depth(self, f):
        '''calculate skin depth'''
        return 1/np.sqrt(MUR0*self.mur*self.sigma*pi*f)

    def periodicity_id(self, b):
        '''identify the periodicity of a given signal'''
        bxy = np.array(b['bxyl'])
        npos = bxy.shape[2]
        if self.symmetry:
            nels = bxy.shape[1]
            ff = []

            for i in range(npos):
                r = npos//(i + 1)

                if r > 1:
                    ff.append(np.sum([np.sum(
                        np.abs(bxy[:, :, :(i+1)]
                               - bxy[:, :, (j+1)*(i+1)-1:(j+2)*(i+1)-1]),
                        axis=1) / ((i+1)*nels)
                                      for j in range(r-1)])/(r-1))
            minf = np.amin(ff)
            i = np.argmin(ff)
            bxymax = np.amax(np.abs(bxy))
            logger.debug("ll %d idx %d minf %f i %d bxymax %d",
                        nels, npos, minf, i, bxymax)

            if minf < bxymax * 1e-4:
                ishift = i
                num_period = (npos - 1) // ishift

                for j in range(num_period - 1):
                    bxy[:, :, 0:i+1] += bxy[:, :, i+j*ishift:(i+1)+(j+1)*ishift]

                bxy[:, :, 0:i + 1] /= num_period
                npos = i + 1
            else:
                ff = [np.sum(np.abs(bxy[:, :, 0] - bxy[:, :, i]))/nels
                      for i in range(npos - 1, 0, -1)]
                minf = np.amin(ff)
                i = np.argmin(ff)
                npos = npos - i

        bx_fft = fft(b['bxyf'][0], (npos-1)/self.tgrid)
        by_fft = fft(b['bxyf'][1], (npos-1)/self.tgrid)

        if self.symmetry:
            bxy_amp = bx_fft.amp + by_fft.amp
            tmp_period = np.array([(npos-1)/i
                                   for i in range((npos-1)//2 + 1) if i > 0])
            npos_nonzero = np.argwhere(tmp_period > 0.1*np.amax(bxy_amp)).squeeze()
            period = tmp_period[npos_nonzero]

            if np.sum(np.around([period[0]%i for i in period])) == 0:
                npos = min(int(np.ceil(np.amax(period))+1), bxy.shape[2])

        try:
            self.tgrid = 60/self.speed*(self.theta[npos-1] - self.theta[0])/360
        except AttributeError:
            self.tgrid = (self.displ[npos-1] - self.displ[0])/self.speed
        logger.debug("Tgrid %f npos %d bx %f by %f xfreq %s yfreq %s",
                     self.tgrid, npos-1,
                     np.max(bx_fft.amp), np.max(by_fft.amp),
                     bx_fft.freq, by_fft.freq)
        return [npos, bx_fft, by_fft]

    def consider_bx(self, wm, hm, bx_fft, by_fft):
        '''check if a caculation is necessary for the x direction'''
        fft_freq = bx_fft.freq
        fft_freq[fft_freq==0] = 0.5e-2

        # skin depth
        delta = self.skin_depth(fft_freq)

        def ampf(bfft, krf):
            return np.sum(bfft.amp**2*bfft.freq*krf)

        if hm > self.lm:
            krfx = dfac(self.lm/delta)
            px = ampf(bx_fft, krfx)*self.lm**3*wm
        else:
            krfx = dfac(hm/delta)
            px = ampf(bx_fft, krfx)*hm**3*wm

        if wm > self.lm:
            krfy = dfac(self.lm/delta)
            py = ampf(by_fft, krfy)*self.lm**3*hm
        else:
            krfy = dfac(wm/delta)
            py = ampf(by_fft, krfy)*wm**3*hm

        if py > 0:
            if px/py > 0.005:
                self.is_x = True

    def bpm_fft(self, nx, ny, nt, elxy, bxy):
        '''interpolate the flux density'''
        x = np.linspace(np.amin(elxy['excpl']),np.amax(elxy['excpl']), nx)
        y = np.linspace(np.amin(elxy['eycpl']),np.amax(elxy['eycpl']), ny)

        xx, yy= np.meshgrid(x, y)
        xx_, yy_ = xx.ravel(), yy.ravel()

        bx_3d = np.zeros((ny,nx,nt))
        by_3d = np.zeros((ny,nx,nt))

        # regular grid
        if self.is_x:
            bx_3d = binterp(elxy['excpl'], elxy['eycpl'],
                                        xx_, yy_, bxy['bxl'][:, 0:nt])
            by_3d = binterp(elxy['excpl'], elxy['eycpl'],
                                        xx_, yy_, bxy['byl'][:, 0:nt])
        else:
            by_3d = binterp(elxy['excpl'], elxy['eycpl'],
                                        xx_, yy_, bxy['byl'][:, 0:nt])

        lfft = (nt-1)//2+1
        by_fft = 2*np.abs(np.fft.rfftn(by_3d[:,:,0:-1]))/(nx*ny*(nt-1))

        py_se = np.zeros((ny,nx,lfft))
        px_se = np.zeros_like(py_se)

        if self.is_x:
            bx_fft = 2*np.abs(np.fft.rfftn(bx_3d[:,:,0:-1]))/(nx*ny*(nt-1))
        else:
            bx_fft = np.zeros_like(by_fft)

        return [bx_fft, by_fft, px_se, py_se]

    def calc_pvpm(self, bamp, f, nu, wm, hm, delta_eff):
        '''calculate eddy current losses for each frequency'''
        wm_nu = wm
        k_xi = 1
        k_eta = 1

        if nu > 1:
            wm_nu = 0.5*wm/(nu - 1)
            k_xi = 0.895
            k_eta = 1.15

        alpha = (hm + delta_eff*self.mur)/hm
        delta = np.sqrt(alpha)*self.skin_depth(f)
        xi = k_xi*wm_nu/delta
        eta = self.lm/(wm_nu*k_eta)

        # caclulation correction factor geometry
        if xi*eta < 500.:
            c_ef_n0 = 32/pi**5/eta/dfac(xi)*6
            sum_r = 0.0
            sum_i = 0.0
            n = np.linspace(100)
            lambda_n = (2*n + 1)*pi
            beta_n = np.sqrt(lambda_n ** 2 + 2j*xi**2)
            beta_nr = np.real(beta_n)
            beta_ni = np.imag(beta_n)

            deno = 1. / ((2*n + 1)**5*abs(beta_n)**6*(cosh(beta_nr*eta) + cos(beta_ni*eta)))
            add_i = ((lambda_n**2 - 2*beta_ni**2)*beta_nr*lambda_n**3*sinh(beta_nr*eta))*deno
            add_r = ((lambda_n**2 + 2*beta_nr**2)*beta_ni*lambda_n**3*sin(beta_ni*eta))*deno

            sum_r = np.sum(np.nan_to_num(add_r, copy=True, nan=0))
            sum_i = np.sum(np.nan_to_num(add_i, copy=True, nan=0))
            c_ef = 1 - c_ef_n0*(sum_i + sum_r)

        else:
            c_ef = 1 + 1/eta

        # calculation correction factor reaction field
        c_rf = dfac(xi)
        if np.isnan(c_rf):
            c_rf = 0

        result = self.sigma*f**2*bamp**2*wm**3*hm*self.lm*c_ef*c_rf
        # calculation ec losses
        if nu > 1:
            return 0.5/(nu-1)**2*result
        else:
            return pi**2/6*result

    def loss(self, bx_fft, by_fft, px_se, py_se, wm, hm):
        '''calculate losses in x and y direction'''
        (ny, nx, lfft) = px_se.shape
        pec = np.zeros((len(self.segx), len(self.segz)))
        nu = np.abs(np.fft.fftfreq(nx, 1/nx))+1
        mu = np.abs(np.fft.fftfreq(ny, 1/ny))+1
        for jj in range(len(self.segx)):
            for kk in range(len(self.segz)):
                for c in range(lfft):
                    for iy in range(ny):
                        for ix in range(nx):
                            if self.is_x and nu[ix] < 2:
                                with warnings.catch_warnings():
                                    warnings.simplefilter('ignore')
                                    px_se[iy,ix,c] = self.calc_pvpm(
                                        bx_fft[iy,ix,c], max(c/self.tgrid, 1e-6),
                                        mu[iy], hm, wm/self.segx[jj], 0)

                            if mu[iy] < 2:
                                with warnings.catch_warnings():
                                    warnings.simplefilter('ignore')
                                    py_se[iy,ix,c] = self.calc_pvpm(
                                        by_fft[iy,ix,c], max(c/self.tgrid, 1e-6),
                                        nu[ix], wm/self.segx[jj], hm, 0)
                py_sum = np.sum(py_se)
                px_sum = np.sum(px_se)
                pec[jj,kk] = (py_sum + px_sum)*(self.ls/self.lm)*self.numpoles*self.segx[jj]

        return np.sum(pec)

    def calc_losses(self):
        '''calculate magnet losses for every load case

        Returns
        --------------
        all_load_cases: list of losses for all load cases
        '''
        all_load_cases = []
        for k in self.pm:
            ialh_loss = 0
            loss_detail = []
            for i in k:
                logger.info('magnet geom / mm: w %.2f h %.2f l %.2f',
                            i["wm"]*1e3, i["hm"]*1e3, i["lm"]*1e3)
                [nt, bx_fft, by_fft] = self.periodicity_id(i['bl'])
                [nx, ny] = ngrid(i['wm'], i['hm'], i['elcp'])
                self.consider_bx(i['wm'], i['hm'], bx_fft, by_fft)
                bfft = self.bpm_fft(nx, ny, nt, i['elcp'], i['bl'])
                loss = self.loss(*bfft, i['wm'], i['hm'])
                ialh_loss += loss
                loss_detail.append([i['spel_key'], loss/self.numpoles])
            self.th_loss.append(loss_detail)
            all_load_cases.append(ialh_loss)

        return all_load_cases



    def ialh2(self, sx_ampl, sy_ampl, sx_phase, sy_phase, freq, wm, hm, lm, delta_eff):
        ''' Calculates eddy current losses for each point (sx, sy, f) in one magnet segment, using ialh2 method

        Inputs:  sx_ampl - dBx/dt amplitudes for each element, for each frequency, for that magnet segment
                 sy_ampl - dBy/dt amplitudes for each element, for each frequency, for that magnet segment
                 sx_phase - x phase of each element for each frequency, for that magnet segment
                 sy_phase - y phase of each element for each frequency, for that magnet segment
                 freq - every considered frequency in Hertz
                 nx, ny - number of elements in each magnet segment
                 wm,hm,lm - width (x), height (y), axial length (z) of each magnet segment
                 delta_eff - effective airgap between magnet and rotor iron
        Returns: Ploss_x - losses [W] generated from x component of magnetic flux, for a given magnet segment
                 Ploss_y - losses [W] generated from y component of magnetic flux, for a given magnet segment
        '''
        sigma_x = self.sigma
        sigma_y = self.sigma
        sigma_z = self.sigma
        delta_xy_x = delta_eff
        delta_xy_y = delta_eff
        delta_yz_y = delta_eff
        delta_xz_x = delta_eff
        epsilon_z_y = 1. + self.mur*delta_xy_y / hm
        epsilon_z_x = 1. + self.mur*delta_xy_x / hm
        epsilon_x_y = 1. + self.mur*delta_yz_y / hm
        epsilon_y_x = 1. + self.mur*delta_xz_x / hm
        rx = sigma_z/sigma_x
        ry = sigma_z/sigma_y
        p0 = 0.5*hm*lm*wm**3*sigma_z

        (nx, ny, nf) = sx_ampl.shape

        Ploss_x = np.zeros((nx,ny,nf))
        Ploss_y = np.zeros_like(Ploss_x)

        for f in range (nf):   #loop over each frequency
            epsilon = wm*np.sqrt(MUR0*self.mur*sigma_z/epsilon_z_y*np.pi*freq[f])

            for n in range (int(nx/2)):
                for m in range (int(ny/2)):
                    symmetry = 4 if (m > 0) & (n > 0) else 2 # symmetry factor 4 due to double symmetry utilization.

                    for k in range (0, 1000, 2): # loop is symmetrical on one side
                        complex_xnmkf = sx_ampl[n,m,f]*2/(np.pi*(k+1))*epsilon_z_x/epsilon_z_y*(m*np.pi*wm/hm)/(2j *epsilon**2 + (n*np.pi)**2 + (m*np.pi*wm/hm)**2 + ry*((k+1)*np.pi*wm/lm)**2)
                        complex_ynmkf = sy_ampl[n,m,f]*2/(np.pi*(k+1)) *n*np.pi/(2j *epsilon**2 + (n*np.pi)**2 + (m*np.pi*wm/hm)**2 + rx*((k+1)*np.pi*wm/lm)**2)
                        real_x_nmkf = np.real(complex_xnmkf)
                        imag_x_nmkf = np.imag(complex_xnmkf)
                        real_y_nmkf = np.real(complex_ynmkf)
                        imag_y_nmkf = np.imag(complex_ynmkf)

                        # Loss component x
                        Plossxnmkf = 2.*symmetry*p0*((sx_ampl[n,m,f]*2/(np.pi*(k+1))*epsilon_y_x/epsilon_z_y)**2*ry
                                                     *((k+1)*np.pi*wm/lm)**2/(4*epsilon**4 + ((n*np.pi)**2 + (m*np.pi*wm/hm)**2 + ry*((k+1)*np.pi*wm/lm)**2)**2) + real_x_nmkf**2 + imag_x_nmkf**2
                                                     - 2*np.cos(sy_phase[n,m,f] - sx_phase[n,m,f])*(real_y_nmkf*real_x_nmkf + imag_y_nmkf*imag_x_nmkf)
                                                     - 2*np.sin(sy_phase[n,m,f] - sx_phase[n,m,f])*(real_y_nmkf*imag_x_nmkf - imag_y_nmkf*real_x_nmkf))
                        # Loss component y
                        Plossynmkf = 2.*symmetry*p0*((sy_ampl[n,m,f]*2/(np.pi*(k+1))*epsilon_x_y/epsilon_z_y)**2*rx
                                                     *((k+1)*np.pi*wm/lm)**2/(4*epsilon**4 + ((n*np.pi)**2 + (m*np.pi*wm/hm)**2 + rx*((k+1)*np.pi*wm/lm)**2)**2) + real_y_nmkf**2 + imag_y_nmkf**2)

                        if (Ploss_x[n,m,f] + Ploss_y[n,m,f]) == 0:      # preventing division by zero in termination criteria evaluation
                            Ploss_x[n,m,f] += Plossxnmkf
                            Ploss_y[n,m,f] += Plossynmkf
                            continue
                        # termination criteria for k loop -> amplitude proportional to 1/k^2
                        if (k > 1) & ((Plossxnmkf + Plossynmkf)/(Ploss_x[n,m,f] + Ploss_y[n,m,f]) < 1e-4):
                            Ploss_x[n,m,f] += Plossxnmkf
                            Ploss_y[n,m,f] += Plossynmkf
                            break

                        Ploss_x[n,m,f] += Plossxnmkf
                        Ploss_y[n,m,f] += Plossynmkf

        return np.sum(Ploss_x), np.sum(Ploss_y)

    def diffDQ(self, bx_pm_3D, by_pm_3D, T):
        ''' Calculates the time derivative of Bx, By

        Inputs:  bx_pm_3D - bx (AC component) for each element for one magnet segment, for each simulation step
                 by_pm_3D - by (AC component) for each element for one magnet segment, for each simulation step
                 T - total simulation time for periodic simulation (see periodicity_id fcn for more)
        Returns: sx_pm - dBx/dt for each element for one magnet segment, for each simulation step
                 sy_pm - dBy/dt for each element for one magnet segment, for each simulation step
        '''
        (nx, ny, nt) = bx_pm_3D.shape
        sx_pm_3D = np.zeros((nx, ny, nt))
        sy_pm_3D = np.zeros_like(sx_pm_3D)

        ti = np.linspace(0, T, nt)
        timestep = ti[1] - ti[0]

        sx_pm_3D = -np.diff(bx_pm_3D, n=1, axis=-1, append=np.reshape(bx_pm_3D[:,:,1], (nx, ny, 1)))/timestep
        sy_pm_3D = -np.diff(by_pm_3D, n=1, axis=-1, append=np.reshape(by_pm_3D[:,:,1], (nx, ny, 1)))/timestep

        return sx_pm_3D, sy_pm_3D

    def diffFreqD(self, bx_pm_3D, by_pm_3D, T):
        ''' Calculates the time derivative of Bx, By, advanced method

        Inputs:  bx_pm_3D - bx (AC component) for each element for one magnet segment, for each simulation step
                 by_pm_3D - by (AC component) for each element for one magnet segment, for each simulation step
                 T - total simulation time for periodic simulation (see periodicity_id fcn for more)
        Returns: sx_pm - dBx/dt for each element for one magnet segment, for each simulation step
                 sy_pm - dBy/dt for each element for one magnet segment, for each simulation step
        '''

        (nx, ny, nt) = bx_pm_3D.shape
        sx_pm_3D = np.zeros((nx, ny, nt))
        sy_pm_3D = np.zeros((nx, ny, nt))
        ti = np.linspace(0, T, nt)
        timestep = ti[1] - ti[0]

        freq = np.fft.rfftfreq(nt-1, timestep)
        nf = freq.shape[0]
        amplbx = np.zeros((freq.shape))
        amplby = np.zeros((freq.shape))
        complbx = np.zeros((nx, ny, nf)).astype(complex)
        complby = np.zeros((nx, ny, nf)).astype(complex)

        for ii in range(nx):
            for jj in range (ny):
                complbx[ii,jj,:] = np.fft.rfftn(bx_pm_3D[ii,jj,0:nt-1])
                complby[ii,jj,:] = np.fft.rfftn(by_pm_3D[ii,jj,0:nt-1])
                amplbx = amplbx + abs(complbx[ii,jj,0:nf])/nf
                amplby = amplby + abs(complby[ii,jj,0:nf])/nf

        amplbx = amplbx/(nx*ny)
        amplby = amplby/(nx*ny)
        amplb = np.sqrt(amplbx**2 + amplby**2)
        fmax2 = 0.5*freq[-1]

        if sum(amplb) == 0:
            warnings.warn('Bx and By data equals to zero, check simulation parameters for this loadcase')
            filt = 0

        if sum(amplb) > 0:
            pec = (np.multiply(amplb,freq))**2
            pecmax = np.max(pec)
            pec = pec/pecmax

            fecmax = fmax2
            imax = int(np.floor(nf/2))
            filt0 = np.ones((nf))
            for ii in range(imax,nf):
                filt0[ii] = fecmax/freq[ii]

            # determine the last significant frequency
            pecf = pec*filt0
            ilim = 0
            feclim = 0
            Ath = 0.05
            for ii in range(1,nf):
                jj = nf - 1 - ii
                if pecf[jj] - pecf[jj + 1] > Ath:
                    ilim = jj
                    feclim = freq[jj]
                    break

            filt = np.ones((nf))
            for ii in range(ilim, nf):
                filt[ii] = feclim/freq[ii]
        for ii in range(nx):       # Derivation in frequency domain
            for jj in range(ny):
                complbx[ii,jj,:] = -complbx[ii,jj,:]*freq*filt*np.pi*2j
                complby[ii,jj,:] = -complby[ii,jj,:]*freq*filt*np.pi*2j

        for ii in range (nx):       # Inverse Fourier-Transformation
            for jj in range (ny):
                sx = np.fft.irfftn(complbx[ii,jj,:], s=[nt - 1], axes=[0])
                sy = np.fft.irfftn(complby[ii,jj,:], s=[nt - 1], axes=[0])
                sx = np.append(sx, sx[0])
                sy = np.append(sy, sy[0])
                sx_pm_3D[ii,jj,:] = sx
                sy_pm_3D[ii,jj,:] = sy

        return sx_pm_3D, sy_pm_3D

    def Process_B_data(self, nx, ny, nsegx, nsegy, nt, elxy, bxy, excpl_new, eycpl_new):
        ''' Processes flux density data: interpolates Bx, By to new resolution defined in Segmentation fcn
                      calculates the dB/dt for x,y axes
                      calculates the FFT of those derivations

        Inputs:  nx,ny - number of elements for each magnet segment (inherited from Segmentation fcn)
                 nsegx,nsegy - number of magnet segments in x,y axis
                 nt - number of time steps (inherited from result of periodicity function), corresponds to 1 period
                 elxy - dictionary with original excpl,eycpl
                 bxy - dictionary with original flux denisities bxl, byl - bx, by in local reference frame
                 excpl_new, eycpl_new - x,y coordinates for new elements (inherited from Segmentation fcn)
        Returns: sx_abs - dBx/dt amplitudes for each element, for each frequency, for entire magnet
                 sy_abs - dBy/dt amplitudes for each element, for each frequency, for entire magnet
                 sx_phase - x phase of each element for each frequency, for entire magnet
                 sy_phase - y phase of each element for each frequency, for entire magnet
                 freq_range - every considered frequency in Hertz
        '''
        nx_tot = int(nx*nsegx)
        ny_tot = int(ny*nsegy)
        bxyl = np.asarray(bxy['bxyl'])
        Bxl_ac = np.zeros((bxyl[0].shape))
        Byl_ac = np.zeros_like(Bxl_ac)

        # Remove the DC component of the original bxl, byl
        for ii in range(Bxl_ac.shape[0]):
            Bxl_ac[ii,:] = bxyl[0, ii, :] - np.mean(bxyl[0,ii,:])
        for ii in range(Byl_ac.shape[0]):
            Byl_ac[ii,:] = bxyl[1,ii,:] - np.mean(bxyl[1,ii,:])
        xx_ = excpl_new.ravel()
        yy_ = eycpl_new.ravel()
        bx_3d_ac = np.zeros((nx_tot,ny_tot,nt))
        by_3d_ac = np.zeros_like(bx_3d_ac)

        ecpl = elxy['ecpl']
        # Interpolation to the new resolution -> [nx*nsegx, ny*nsegy, nt]
        by_3d_ac = binterp_ialh2(ecpl[0], ecpl[1],
                                 xx_, yy_, Byl_ac[:, 0:nt])
        if self.is_x:
            bx_3d_ac = binterp_ialh2(ecpl[0], ecpl[1],
                                     xx_, yy_, Bxl_ac[:, 0:nt])
        bx_3d_ac = bx_3d_ac.reshape(nx_tot,ny_tot,nt)
        by_3d_ac = by_3d_ac.reshape(nx_tot,ny_tot,nt)

        xx_ = xx_.reshape(nx_tot, ny_tot)
        yy_ = yy_.reshape(nx_tot, ny_tot)

        if nt % 2:
            nf = int((nt - 1)/2 + 1)
        else:
            nf = int(nt / 2)

        sx_abs = np.zeros((2*nx_tot, 2*ny_tot, nf))
        sy_abs = np.zeros_like(sx_abs)
        sx_phase = np.zeros_like(sx_abs)
        sy_phase = np.zeros_like(sx_phase)
        for ii in range (nsegx):
            for jj in range (nsegy):

                diffgrad = 3   # choose the derivation metod. diffFreqD method recommended
                if diffgrad == 1:
                    sx_pm, sy_pm = self.diffDQ(bx_3d_ac[ii*nx:(ii+1)*nx, jj*ny:(jj+1)*ny,:],
                                               by_3d_ac[ii*nx:(ii+1)*nx, jj*ny:(jj+1)*ny,:],
                                               self.tgrid)
                else:
                    sx_pm, sy_pm = self.diffFreqD(bx_3d_ac[ii*nx:(ii+1)*nx, jj*ny:(jj+1)*ny,:],
                                                  by_3d_ac[ii*nx:(ii+1)*nx, jj*ny:(jj+1)*ny,:],
                                                  self.tgrid)

                Gx_seg = np.zeros((2*nx, 2*ny, nt-1))                  # omit the last step of the time vector -> duplicated with the first step
                Gx_seg[0:nx, 0:ny,:] = +sx_pm[:,:, 0:-1]
                Gx_seg[0:nx,ny:,:] = -sx_pm[:,::-1, 0:-1]          # this section flips and "doubles" the data for each segment, so FFT can include every point properly
                Gx_seg[nx:,:,:] = np.flipud(Gx_seg[0:nx,:,:])

                Gy_seg = np.zeros((2*nx, 2*ny, nt-1))                  # omit the last step of time vector -> duplicated with the first step
                Gy_seg[0:nx, 0:ny,:] = +sy_pm[:,:, 0:-1]
                Gy_seg[nx:,0:ny,:] = -sy_pm[::-1,:, 0:-1]          # this section flips and "doubles" the data for each segment, so FFT can include every point properly
                Gy_seg[:,ny:,:] = np.fliplr(Gy_seg[:,0:ny,:])

                sx_FFT_seg = np.fft.rfftn(Gx_seg)
                sy_FFT_seg = np.fft.rfftn(Gy_seg)
                sx_abs_seg = 2*abs(sx_FFT_seg /(2*nx*2*ny*(nt - 1)))    # dBx/dt amplitudes for each segment element, for each frequency for that magnet segment
                sy_abs_seg = 2*abs(sy_FFT_seg /(2*nx*2*ny*(nt - 1)))    # dBy/dt amplitudes for each segment element, for each frequency for that magnet segment
                sx_phase_seg = np.angle(sx_FFT_seg)
                sy_phase_seg = np.angle(sy_FFT_seg)
                freq_range = np.fft.rfftfreq(nt-1, self.tgrid/(nt-1))   # every considered frequency in Hertz

                sx_abs[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:] = sx_abs_seg     # dBx/dt amplitudes for each element, for each frequency, for entire magnet
                sy_abs[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:] = sy_abs_seg     # dBy/dt amplitudes for each element, for each frequency, for entire magnet
                sx_phase[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:] = sx_phase_seg # x phase of each element for each frequency, for entire magnet
                sy_phase[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:] = sy_phase_seg # y phase of each element for each frequency, for entire magnet

        return sx_abs, sy_abs, sx_phase, sy_phase, freq_range

    def loss_ialh2(self, sx_abs, sy_abs, sx_phase, sy_phase, freq_range, nx, ny, wm, hm, lm, nsegx, nsegy, nsegz, delta_eff):
        ''' Loops over each magnet segment and calculates the losses for the entire magnet

        Inputs:  sx_abs, sy_abs - dBx/dt, dBy/dt amplitudes for each element, for each frequency, in entire magnet
                 sx_phase, sy_phase - corresponding phases of each freq for each element in entire magnet
                 freq_range - every considered frequency in Hz
                 nx, ny - number of elements in each magnet segment in x and y direction
                 wm, hm, lm - total magnet width (x), height (y), axial length (z)
                 nsegx, nsegy, nsegz - number of magnet segments in x,y,z direction
                 delta_eff - needed for ialh2 losses calculation, effective airgap between magnet and rotor iron
        Returns: total eddy current losses for entire magnet (all segments), for 1 magnet of the machine
        '''
        pec = np.zeros((nsegx, nsegy, nsegz))

        for ii in range(nsegx):
            for jj in range(nsegy):    # nsegy is always = 1
                for kk in range(nsegz):
                    Plossx, Plossy = self.ialh2(sx_abs[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:],
                                                sy_abs[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:],
                                                sx_phase[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:],
                                                sy_phase[2*ii*nx:2*(ii+1)*nx, 2*jj*ny:2*(jj+1)*ny,:],
                                                freq_range, wm/nsegx, hm/nsegy, lm/nsegz, delta_eff)
                    pec[ii,jj,kk] = (Plossx + Plossy)

        return np.sum(pec)  # total eddy current losses for the entire magnet (all segments)


    def calc_losses_ialh2(self, nsegx=1, nsegy=1, nsegz=1):
        ''' Calculates magnet losses for every load case

        Inputs:  number of magnet segments in x,y,z direction
        Returns: all_load_cases: list of losses for all load cases
        '''

        nsegx = max(1, nsegx)    #  1 = no segmentation
        nsegz = max(1, nsegz)    #  1 = no segmentation
        nsegy = 1               # y segmentation not supported, nsegy is always = 1

        delta_eff = 0

        all_load_cases = []
        for k in self.pm:               # loop for each load case
            ialh_loss = 0
            loss_detail = []
            for i in k:                 # loop for each superelement in a case
                logger.info('magnet geom / mm: w %.2f h %.2f l %.2f segments %s',
                            i["wm"]*1e3, i["hm"]*1e3, i["lm"]*1e3,
                            (nsegx, nsegy, nsegz))

                (nt, bx_fft, by_fft) = self.periodicity_id(i['bl'])  # finds the time periodic part of the simulation
                (nx, ny, nz, excpl_new, eycpl_new) = Segmentation(
                    i['wm'], i['hm'], i['lm'], i['elcp'],
                    nsegx, nsegy, nsegz)

                wm = i['wm']
                hm = i['hm']
                lm = i['lm']
                self.consider_bx(wm, hm, bx_fft, by_fft)
                (sx_abs, sy_abs, sx_phase, sy_phase, freq_range) = self.Process_B_data(
                    nx, ny, nsegx, nsegy, nt, i['elcp'], i['bl'], excpl_new, eycpl_new)
                loss = self.loss_ialh2(sx_abs, sy_abs, sx_phase, sy_phase, freq_range,
                                       nx, ny, wm, hm, lm, nsegx, nsegy, nsegz, delta_eff) * self.numpoles
                ialh_loss += loss
                logger.info('Loadcase %d, Superelement %s, Total losses =  %.3f W',
                            i["loadcase"], i["spel_key"], loss)
                loss_detail.append([i['spel_key'], loss/self.numpoles])
            self.th_loss.append(loss_detail)
            all_load_cases.append(ialh_loss)

        return all_load_cases
