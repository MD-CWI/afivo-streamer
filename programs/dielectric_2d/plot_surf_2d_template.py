#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 15:34:12 2022

@author: xiaoran
"""


import numpy as np
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

font = {
    "font.family": "serif",
    "font.serif": "Times New Roman",
    "font.size": 14
}

font_cn = {
    "family": "SimSun",
    "size": 14
}
mpl.rcParams.update(font)
plt.rcParams['mathtext.fontset'] = 'cm'  # 'cm' (Computer Modern)

width_savgol_filter = 15
order_savgol_filter = 1

fname = ['test_forward_euler_000020_data', 'test_heuns_method_000020_data','test_midpoint_method_000020_data']

fig = plt.figure(figsize=(5, 4))
time = np.array([28])
styl_list=['-','--','-.',':','-','--','-.',':']
for fn, styl in zip(fname, styl_list):
    data = np.load(fn+'.npy')
    data = data[np.argsort(data[:,1])]
    data[:,2:] = data[:,2:]*1e6
    data[:,:2] = data[:,:2]*1e3
    
    total = data[:,4]
    total = savgol_filter(total, width_savgol_filter, order_savgol_filter)
    plt.plot(data[:,1], total, label=fn, ls=styl)

# plt.axis([40, 0, -400, 20])
plt.xticks()
plt.yticks()

plt.rcParams['xtick.direction'] = 'in'  # 将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'  # 将y轴的刻度方向设置向内
plt.minorticks_on()
plt.tick_params(top=True, bottom=True, left=True, right=True)
plt.xlabel('y (mm)')
plt.ylabel('$\\sigma_s$ (pC/mm$^2$)')
plt.grid(linestyle='-.', alpha=0.7)
plt.legend(fontsize=14)
#plt.subplots_adjust(left = 0.3,bottom=0.3,right=0.8, top=0.8)
# plt.savefig(fname+'.png',dpi=300,bbox_inches='tight')

# plt.savefig('negative_total_charge_md.png', dpi=300, bbox_inches='tight')