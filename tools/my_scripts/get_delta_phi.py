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
    delta_phi = np.array([])
    phi_0 = np.array([])
    phi_t = np.array([])
    try:
        log = np.loadtxt(f"{f}_log.txt", skiprows=1)
        log_phi0 = np.loadtxt(f"{f}_phi_0000.curve", skiprows=1)
    except IOError:
        break
    log_phi0 = truncate_axis(log_phi0, zmax=args.zmax, col=0)
    z0 = log_phi0[:,0]
    phi0 = log_phi0[:,1]
    f_phi0 = interpolate.UnivariateSpline(z0, phi0, s=0)
    j = args.i0
    while (j >= 0):
        if (j > args.i1):
            break
        try:
            log_Ech = np.loadtxt(f"{f}_electric_fld_{str(j).zfill(4)}.curve", skiprows=1)
            log_phi = np.loadtxt(f"{f}_phi_{str(j).zfill(4)}.curve", skiprows=1)
        except IOError:
            break
        log_Ech = truncate_axis(log_Ech, zmax=args.zmax, col=0)
        z = log_Ech[:,0]
        Ech = log_Ech[:,1]
        k = np.argmax(Ech)
        if (j == 0):
            z_head = log[j,9]
        else:
            z_head = z[k]
        log_phi = truncate_axis(log_phi, zmax=args.zmax, col=0)
        z = log_phi[:,0]
        phi = log_phi[:,1]
        f_phi = interpolate.UnivariateSpline(z, phi, s=0)
        t = np.append(t, log[j,1])
        y =np.append(y, z_head)
        delta_phi = np.append(delta_phi, f_phi(z_head)-f_phi0(z_head))
        phi_0 = np.append(phi_0, f_phi0(z_head))
        phi_t = np.append(phi_t, f_phi(z_head))
        print(f"Processing {f}_phi_{str(j).zfill(4)}.curve")
        j += 1
    np.savetxt(f"y_delta_phi_{f}.txt", np.column_stack([t, y, delta_phi, phi_0, phi_t]), 
               header='time  y  delta_phi  phi_0  phi_t', comments='', delimiter='  ', fmt='%.8e')    
    print(f"Saved y_delta_phi_{f}.txt")
