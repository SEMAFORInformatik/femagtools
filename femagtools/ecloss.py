'''
    femagtools.ecloss
    ~~~~~~~~~~~~~~~~

    Calculate Magnet Losses with IALH Method
'''
__author__ = 'Max Hullmann, Dapu Zhang'

import logging
import warnings
from .amela import Amela
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
    be = np.sqrt(wm*hm/elxy['excp'].shape[0])
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
    f = RBFInterpolator(np.array([[i, j] for i, j in zip(x, y)]), b) 
    inp = f(np.array([[i, j] for i, j in zip(xq, yq)]))
    return inp.reshape(len(np.unique(yq)), len(np.unique(xq)))


class MagnLoss(Amela): 
    '''Calculate Magnet Losses with IALH Methode
    Parameters
    ----------
    workdir: working directory 
    modelname: name of the femag model (*.nc)
    ibeta: load cases [0, 1, 2]
    '''
    def __init__(self, workdir, modelname, ibeta, **kwargs): 
        super().__init__(workdir, magnet_data=dict(name=modelname))
        self.pm = self.get_magnet_data_all(ibeta)
        self.theta = self.pm[-1][-1]['phi'] # rotor pos
        self.speed = kwargs.get('speed', self.pm[-1][-1]['speed'])
        self.tgrid = 60/self.speed*(self.theta[-1] - self.theta[0])/360
        self.lt = len(self.theta)
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
        self.is_meter = False
        # determine the number of segments in z direction
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
        bx = b['bxl']
        by = b['byl']
        idx = bx.shape[1]
        if self.symmetry:
            ll = bx.shape[0]
            ff = []

            for i in range(idx):
                r = idx//(i + 1)

                if r > 1:
                    f = 0

                    for j in range(r-1):
                        f += np.sum(np.abs(bx[:, 0:(i+1)] - bx[:, (j+1)*(i+1)-1:(j+2)*(i+1)-1]))/((i+1)*ll)
                        f += np.sum(np.abs(by[:, 0:(i+1)] - by[:, (j+1)*(i+1)-1:(j+2)*(i+1)-1]))/((i+1)*ll)

                    ff.append(f/(r-1))

            minf = np.amin(ff)
            i = np.argmin(ff)
            bxymax = np.amax([np.amax(np.abs(bx)), np.amax(np.abs(by))])

            if minf < bxymax * 1e-4:
                ishift = i
                num_period = (idx - 1) // ishift

                for j in range(num_period - 1):
                    bx[:, 0:i+1] += bx[:, i+j*ishift:(i+1)+(j+1)*ishift]
                    by[:, 0:i+1] += by[:, i+j*ishift:(i+1)+(j+1)*ishift]

                bx[:, 0:i + 1] /= num_period
                by[:, 0:i + 1] /= num_period
                idx = i + 1
            else:
                ff = []

                for i in range(idx - 1, 0, -1):
                    f1 = np.sum(np.abs(bx[:, 0] - bx[:, i]))/ll + np.sum(np.abs(by[:, 0] - by[:, i]))/ll
                    ff.append(f1)

                minf = np.amin(ff)
                i = np.argmin(ff)
                idx = idx - i
    
        bx_fft = fft(b['bxf'][0:idx-1], (idx-1)/self.tgrid)
        by_fft = fft(b['byf'][0:idx-1], (idx-1)/self.tgrid)
        
        if self.symmetry: 
            bxy_amp = bx_fft.amp + by_fft.amp
            tmp_period = np.array([(idx-1)/i for i in range((idx-1)//2 + 1) if i > 0])
            idx_nonzero = np.argwhere(tmp_period > 0.1*np.amax(bxy_amp)).squeeze()
            period = tmp_period[idx_nonzero]

            if np.sum(np.around([period[0]%i for i in period])) == 0: 
                idx = int(np.ceil(np.amax(period))+1)
                if idx > bx.shape[1]: 
                    idx = bx.shape[1]  

        self.tgrid = 60/self.speed*(self.theta[idx-1] - self.theta[0])/360

        return [idx, bx_fft, by_fft]

    def consider_bx(self, wm, hm, bx_fft, by_fft):
        '''check if a caculation is necessary for the x direction'''
        fft_freq = bx_fft.freq
        fft_freq[fft_freq==0] = 0.5e-2

        if not self.is_meter: 
            self.ls *= 1e-3
            self.lm *= 1e-3
            self.is_meter = True
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

        if px/py > 0.005:
            self.is_x = True
        return ' '
        
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
            for t in range(nt): 
                bx_3d[:, :, t] = binterp(elxy['excpl'], elxy['eycpl'], 
                                         xx_, yy_, bxy['bxl'][:, t])
                by_3d[:, :, t] = binterp(elxy['excpl'], elxy['eycpl'], 
                                         xx_, yy_, bxy['byl'][:, t])
        else: 
            for t in range(nt): 
                by_3d[:, :, t] = binterp(elxy['excpl'], elxy['eycpl'], 
                                         xx_, yy_, bxy['byl'][:, t])
             
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
            n = np.array([i for i in range(100)])
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
                                    px_se[iy,ix,c] = self.calc_pvpm(bx_fft[iy,ix,c], max(c/self.tgrid, 1e-6),
                                                                    mu[iy], hm, wm/self.segx[jj], 0)
                                
                            if mu[iy] < 2:
                                with warnings.catch_warnings():
                                    warnings.simplefilter('ignore')
                                    py_se[iy,ix,c] = self.calc_pvpm(by_fft[iy,ix,c], max(c/self.tgrid, 1e-6),
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
                logger.info(f'magnet width and height: {i["wm"]:.2f}mm {i["hm"]:.2f}mm')
                [nt, bx_fft, by_fft] = self.periodicity_id(i['bl'])
                [nx, ny] = ngrid(i['wm'], i['hm'], i['elcp'])
                keyset = ('wm', 'hm')
                for j in keyset: 
                    i[j]*=1e-3
                self.consider_bx(i['wm'], i['hm'], bx_fft, by_fft)
                bfft = self.bpm_fft(nx, ny, nt, i['elcp'], i['bl'])
                loss = self.loss(*bfft, i['wm'], i['hm'])
                ialh_loss += loss
                loss_detail.append([i['spel_key'], loss/self.numpoles])
            self.th_loss.append(loss_detail)
            all_load_cases.append(ialh_loss)

        return all_load_cases
