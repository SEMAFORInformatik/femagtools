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


def harm_leak_ind(E_fft, order_fft, freq, Ia): # needs to be validated
    '''returns harmonic leakage inductance per phase'''
    L_harm = []
    for ii in range(1, len(order_fft)): # skip 1st order
        L_harm.append(E_fft[ii] / Ia / (2*np.pi*freq*order_fft[ii]))

    return sum(L_harm)

