# -*- coding: utf-8 -*-
'''Calculate leakage inductances'''
import numpy as np
import logging

logger = logging.getLogger(__name__)


def end_wdg_leak_ind_round_wires(l_ew, p, q, m, y, taup, w1, layers, num_par_wdgs, slot_height):
    '''returns end winding leakage inductance per phase per end winding side (DE or NDE)'''
    mue0 = 4*np.pi*1e-7
    if layers == 2:
        lambda_ew = 0.34*q*(1 - 2*y*taup/(np.pi*l_ew*m*q))
    if layers == 1:
        lambda_ew = q*(0.67 - 0.43*(taup + slot_height*np.pi/(2*p))/l_ew)
    L_ew = mue0*2/(p*q)*(w1/num_par_wdgs)**2*l_ew*lambda_ew
    return L_ew


def end_wdg_leak_ind_hairpin_wires(): #TODO
    L_ew = 0
    return L_ew


#def slot_leakage_inductance_round_wires(p, q, w1, num_par_wdgs, layers):
#    '''returns slot leakage inductance per phase'''
#    mue0 = 4*np.pi*1e-7
#    if layers == 1:
#        lambda_slot = 0 # tbc
#    if layers == 2:
#        t1 = b2/bs
#        t2 = b1/b2
#        t12 = b1/b2
#        kt1 = (4*t1**2 - t1**4 - 4*np.log(t1) -3)/(4*(1 - t1)*(1 - t1**2)**2) if t1 != 1 else 0
#        kt2 = (4*t2**2 - t2**4 - 4*np.log(t2) - 3)/(4*(1 - t2)*(1 - t2**2)**2) if t2 != 1 else 0
#        kt12 = (t12**2 - 2*np.log(t12) - 1)/(2*(1 - t12)*(1 - t12**2)) if t12 != 1 else 0
#        const = 0.1424 + 0.5*np.arcsin(np.sqrt(1 - (bo/b1)**2)) + ho/bo
#        lambda_t = h2/b2*kt2 + const
#        lambda_b = h3/bs*kt1 + h2/(b2-b1)*np.log(b2/b1) + const if b2 != b1 else h3/bs*kt1 + const
#        lambda_tb = h2/b2*kt12 +  const
#        lambda_slot = lambda_t + lambda_b + lambda_tb
#    L_slot = mue0*2/(p*q)*(w1/num_par_wdgs)**2*lambda_slot
#    return L_slot


def slot_leak_ind_fea(): #TODO
    '''returns slot and tooth tip leakage inductance'''
    # make a single slot model with detailed windings
    # run current through windings
    # L_slot = flux / current
    # scale to get values per phase
    L_slot = 0
    return L_slot


def harm_leak_ind(E_fft, order_fft, freq, Ia): # needs to be validated
    '''returns harmonic leakage inductance per phase'''
    L_harm = []
    for ii in range(1, len(order_fft)): # skip 1st order
        L_harm.append(E_fft[ii] / Ia / (2*np.pi*freq*order_fft[ii]))

    return sum(L_harm)

