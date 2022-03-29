#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:32:05 2021

@author: Baohong.Guo
"""

import numpy as np
import argparse
from scipy import interpolate

def get_descend_Emax_max_index(Emax, Emax_min=3.5e6):  
    index = np.size(Emax)-1
    for i in range(np.size(Emax)):
        if (Emax[i] < Emax_min):
            index = i-1
            break
    return index   

def truncate_axis(array, zmax=0.1151, zmin=0, col=0): 
    array = array.astype(float)
    z = array[:,int(col)]
    i1 = 0
    i2 = np.size(z)
    for i in range(np.size(z)):
        if z[i] >= zmin:
            i1 = i
            break
    for i in range(np.size(z)):
        if z[i] > zmax:
            i2 = i
            break
    return array[i1:i2]

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("filename", type=str, nargs='+', help="Input filename(s)")
p.add_argument('-zmax', type=float, default=0.1151, help='z value at the electrode tip')
args = p.parse_args()

filename = np.array([])
it = np.array([]).astype(int)
t = np.array([])
y = np.array([])
Est_dV_t = np.array([])
Est_dV_0 = np.array([])
L = np.array([])
delta_phi_t = np.array([])
delta_phi_0 = np.array([])
phi_head_t = np.array([])
phi_head_0 = np.array([])
phi_tail = np.array([])
y_tail = np.array([])
    
for i, f in enumerate(args.filename):
    try:
        log = np.loadtxt(f"{f}_log.txt", skiprows=1)
        Emax = log[:,7]
        j = get_descend_Emax_max_index(Emax, Emax_min=3.5e6) 
        log_Ebg = np.loadtxt(f"{f}_electric_fld_0000.curve", skiprows=1)
        log_phi0 = np.loadtxt(f"{f}_phi_0000.curve", skiprows=1)
        log_phi = np.loadtxt(f"{f}_phi_{str(j).zfill(4)}.curve", skiprows=1)
    except IOError:
        break
    log_Ebg = truncate_axis(log_Ebg, zmax=args.zmax, col=0)
    zbg = log_Ebg[:,0]
    Ebg = log_Ebg[:,1]
    k = np.argmax(Ebg)
    y_end = zbg[k] 
    y_head = log[j,9]
    Ls = np.abs(y_head-y_end)
    log_phi0 = truncate_axis(log_phi0, zmax=args.zmax, col=0)
    z0 = log_phi0[:,0]
    phi0 = log_phi0[:,1]
    f_phi0 = interpolate.UnivariateSpline(z0, phi0, s=0)
    log_phi = truncate_axis(log_phi, zmax=args.zmax, col=0)
    z = log_phi[:,0]
    phi = log_phi[:,1]
    f_phi = interpolate.UnivariateSpline(z, phi, s=0)
    
    filename = np.append(filename, f)
    it = np.append(it, int(log[j,0]))
    t = np.append(t, log[j,1])
    y =np.append(y, y_head)
    Est_dV_t = np.append(Est_dV_t, np.abs(f_phi(y_head)-f_phi(y_end))/Ls)
    Est_dV_0 = np.append(Est_dV_0, np.abs(f_phi0(y_head)-f_phi0(y_end))/Ls)
    L = np.append(L, Ls)
    delta_phi_t = np.append(delta_phi_t, f_phi(y_head)-f_phi(y_end))
    delta_phi_0 = np.append(delta_phi_0, f_phi0(y_head)-f_phi0(y_end))
    phi_head_t = np.append(phi_head_t, f_phi(y_head))
    phi_head_0 = np.append(phi_head_0, f_phi0(y_head))
    phi_tail = np.append(phi_tail, f_phi0(y_end))
    y_tail = np.append(y_tail, y_end)
    print(f"Processing {f}_phi_{str(j).zfill(4)}.curve")

savename = f[::-1].split("_", 1)[1][::-1]
np.savetxt(f"y_stability_field_{savename}_{i+1}cases.txt", np.column_stack([filename, it, t, y, 
           Est_dV_t, Est_dV_0, L, delta_phi_t, delta_phi_0, phi_head_t, phi_head_0, phi_tail, y_tail]), 
           header='filename  it  t  y  Est_dV_t  Est_dV_0  L  delta_phi_t  delta_phi_0  phi_head_t  phi_head_0  phi_tail  y_tail', 
           comments='', delimiter='  ', fmt='%s')    
print(f"Saved y_stability_field_{savename}_{i+1}cases.txt")
