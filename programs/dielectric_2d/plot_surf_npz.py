#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 17:38:31 2022

@author: xiaoran
"""

import numpy as np
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

qe = 1.6022e-19

fname = 'test_midpoint_method_000020_surface.npz'
npzfile = np.load(fname)
nc = npzfile['info'][3]
direction = int(npzfile['000001_loc'][-1])
dim = int((np.size(npzfile['000001_loc'])-1)/2)
name_files = npzfile.files
n = int((len(name_files)-1)/2)

data = []
# form a list
# loc(:) density(:)
for i in range(n):
    i = i+1
    loc = str(i).zfill(6)+'_loc'
    sd = str(i).zfill(6)+'_sd'
    for j in range(nc):
        xyz = np.array([npzfile[loc][0], npzfile[loc][1]+j*npzfile[loc][3]])
        surfcharge = np.array(npzfile[sd][j])*qe
        data.append(np.hstack([xyz, surfcharge]))

data = np.array(data)

sname = fname[:-11]+'data'
np.save(sname, data)
print('save:', sname)
