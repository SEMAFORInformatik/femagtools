import numpy as np 
from .utils import fft
import logging 

logger = logging.getLogger('femagtools.semi_fea')

def shift_array(v, idx): 
    '''shift array by index'''
    return v[idx::] + v[0:idx]

def fft_filter(result_fft, perc=0.01):
    '''filter FFT result with amplitude'''
    result = {"order": [], "y":[]}
    base_amp = result_fft['a']
    for i, j in enumerate(result_fft['nue']):
        if j >= perc*base_amp: 
            result['order'].append(i)
            result['y'].append(j)
    return result

def fast_skew_cogg(result, skew_setup): 
    '''Calculate cogging torque/Back-EMF with step skewing based on unskewed result
    Arguments: 
    result: BCH objects 
    skew_setup(dict): {"skew_angle": 10, "nu_skew_steps": 2}
    '''
    skew_angle = skew_setup['skew_angle']
    num_skew_steps =skew_setup['num_skew_steps']
    skew_angle_intern = 0.0
    skew_angle_array = []
    T_slice = []
    bemf = {"1": [], "2": [], "3": []}
    bemf_slice = {"1": [], "2": [], "3": []}
    bemf_skew = {"1": [], "2": [], "3": []}
    bemf_skew_fft = {"1": [], "2": [], "3": []}

    keyset = ('1', '2', '3')

    # check if skew steps equals 2
    if num_skew_steps == 2: 
        skew_angle_intern = skew_angle
        skew_angle_array = [-skew_angle_intern/2, skew_angle_intern/2]
    else: 
        skew_angle_intern = skew_angle/num_skew_steps*(num_skew_steps-1)/2
        skew_angle_array = np.linspace(-skew_angle_intern, skew_angle_intern, 
                                       num_skew_steps, endpoint=True).tolist()

    angle = result.torque[-1]['angle']
    T = result.torque[-1]['torque'][0:-1]
    # get back-emf from BCH
    for i in keyset: 
        bemf[i] = result.flux[i][-1]['voltage_dpsi'][0:-1]

    angl_resl = angle[1]
    tmp = np.unique(np.abs(skew_angle_array))
    skew_angl_resl = 0.0
    if np.amin(tmp) == 0.0: 
        skew_angl_resl = tmp[1]
    else: 
        skew_angl_resl = tmp[0]
    
    divider = skew_angl_resl/angl_resl 
    if divider - np.floor(divider) > 1e-15: 
        # TODO: Interpolation if angle resolution doesn't match
        logger.warning("Wrong Mesh Size in the airgap mesh")
    else: 
        logger.info(f"number of element shifted {divider}")

    for i in skew_angle_array: 
        idx = int(i/angl_resl)
        if i != 0: 
            T_slice.append(shift_array(T, idx))
            for j in keyset:
                bemf_slice[j].append(shift_array(bemf[j], idx))
        else: 
            # do nothing
            T_slice.append(T)
            for j in keyset:
                bemf_slice[j].append(bemf[j])

    # average torque
    T_sum = 0
    for i in T_slice: 
        T_sum += np.array(i)
    T_skew = (T_sum/num_skew_steps).tolist()
    T_skew += [T_skew[0]]
    T_fft = fft_filter(fft(angle, T_skew, pmod=2))

    # average back-emf
    for j in keyset:
        flx_skew = 0
        for k in bemf_slice[j]:
            flx_skew+=np.array(k)
        bemf_skew[j] = (flx_skew/num_skew_steps).tolist()
        bemf_skew[j] += [bemf_skew[j][0]]
        bemf_skew_fft[j] = fft_filter(fft(angle, bemf_skew[j], pmod=2))

    for i in range(len(T_slice)):
        T_slice[i] = (np.array(T_slice[i])/num_skew_steps).tolist()
        T_slice[i]+=[T_slice[i][0]]
    
    return {"angle": angle,
            "cogging_torque": T_skew, 
            "cogging_torque_fft": T_fft, 
            "BEMF": bemf_skew, 
            "BEMF_fft": bemf_skew_fft, 
            "cogging_torque_slice": T_slice}
    