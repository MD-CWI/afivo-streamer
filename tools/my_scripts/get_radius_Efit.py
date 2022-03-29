#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:32:05 2021

@author: Baohong.Guo
"""

import numpy as np
import argparse
from scipy import interpolate

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

def get_radius_Efit(log_Ebg, log_Ech):
    z0 = log_Ebg[:,0]
    Ebg = log_Ebg[:,1]
    f_Ebg = interpolate.UnivariateSpline(z0, Ebg, s=0)
    z = log_Ech[:,0]
    Ech = log_Ech[:,1]
    i = np.argmax(Ech)
    E_R = np.amax(Ech)/4 + 3*f_Ebg(z[i])/4
    E_head = Ech[0:i+1] 
    i0 = np.amin(np.where(E_head >= (0.01*np.amax(E_head)+0.99*np.amin(E_head))))
    i1 = np.amax(np.where(E_head <= (0.99*np.amax(E_head)+0.01*np.amin(E_head))))
    E_head = E_head[i0:i1+1]
    z_head = z[i0:i1+1]
    f_z = interpolate.UnivariateSpline(E_head, z_head, s=0)
    R = z[i]-f_z(E_R)
    return i, z[i], R

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("filename", type=str, nargs='+', help="Input filename(s)")
p.add_argument('-zmax', type=float, default=0.1151, help='z value at the electrode tip')
p.add_argument('-i0', type=int, default=0, help='Start index in database')
p.add_argument('-i1', type=int, default=10000, help='Stop index in database')
args = p.parse_args()

for i, f in enumerate(args.filename):
    t = np.array([])
    y = np.array([])
    R_Efit = np.array([])
    R_Er = np.array([])
    try:
        log = np.loadtxt(f"{f}_log.txt", skiprows=1)
        log_Ebg = np.loadtxt(f"{f}_electric_fld_0000.curve", skiprows=1)
    except IOError:
        break
    log_Ebg = truncate_axis(log_Ebg, zmax=args.zmax, col=0)
    j = args.i0
    while (j >= 0):
        if (j > args.i1):
            break
        try:
            log_Ech = np.loadtxt(f"{f}_electric_fld_{str(j).zfill(4)}.curve", skiprows=1)
        except IOError:
            break
        log_Ech = truncate_axis(log_Ech, zmax=args.zmax, col=0)
        t = np.append(t, log[j,1])
        R_Er = np.append(R_Er, log[j,14])
        if (j == 0):
            y = np.append(y, log[j,9])
            R_Efit = np.append(R_Efit, 0)
        else:
            y = np.append(y, get_radius_Efit(log_Ebg, log_Ech)[1])
            R_Efit = np.append(R_Efit, get_radius_Efit(log_Ebg, log_Ech)[2])
        print(f"Processing {f}_electric_fld_{str(j).zfill(4)}.curve")
        j += 1
    np.savetxt(f"y_radius_Efit_{f}.txt", np.column_stack([t, y, R_Efit, R_Er]), 
               header='time  y  R_Efit  R_Er', comments='', delimiter='  ', fmt='%.8e')    
    print(f"Saved y_radius_Efit_{f}.txt")
