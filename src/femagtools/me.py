"""
Postprocessing FEMAG-ME results
"""
__author__ = 'dapu zhang'

import re 
import numpy as np 
import matplotlib.pyplot as plt
import logging

#  set logging
logger = logging.getLogger(__name__)

def get_eigenvectors(filelist, psfiles):
    """return eigenfrequency and eigenvectors (raw data and psfiles)"""
    eigenvecs = {}
    eigenfreq = []
    ps_data = []
    temp = [] 
    temp_arr = []
    read_fig = False

    if len(psfiles) > 0: read_fig = True
    if len(filelist) == 0: return []
    
    for k, kk in zip(filelist, range(len(filelist))): 
        # read ps file
        if read_fig: ps_data.append(plt.imread(psfiles[kk]))
        eigenvecs[str(kk)] = {}
        with open(k, 'r') as f: 
            data = f.read().split('\n')
        eigenval = float(re.split('\s+', data[0].strip())[2])
        eigenfreq.append(eigenval/2/np.pi)
        
        for i in range(len(data[1:-1])): 
            for j in data[2+i].strip().split(' '): 
                if j != '':
                    temp.append(float(j))

        temp_arr = np.array(temp).reshape(-1, 5)
        idx_collect = []
        # remove the index, at which the dof equals zero
        for c in range(temp_arr.shape[0]): 
            if np.abs(temp_arr[c, 1]) < 1e-15 and \
            np.abs(temp_arr[c, 2]) < 1e-15: 
                idx_collect.append(c)
        # prepare the data from the further calculation
        temp_arr = np.delete(temp_arr, np.array(idx_collect), axis=0)
        temp_arr = np.insert(temp_arr, 3, np.zeros((temp_arr.shape[0], )), axis=1)
        temp_arr = np.hstack((temp_arr, np.zeros((temp_arr.shape[0], 1))))
        eigenvecs[str(kk)].update({'eigenvecs': temp_arr[:, 1:4]})
        # return node position
        if kk == 0: eigenvecs.update({'node_pos': temp_arr[:, 4::]})

    return [ps_data, eigenfreq, eigenvecs] 