#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:32:05 2021

@author: Baohong.Guo
"""

import numpy as np
import argparse
from scipy import integrate

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

def get_q_integral(cross):
    z = cross[:,0]
    λ = cross[:,2]*1.6e-19
    q = integrate.trapz(λ, z)
    return q

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("filename", type=str, nargs='+', help="Input filename(s)")
p.add_argument('-zmax', type=float, default=0.1151, help='zmax')
p.add_argument('-i0', type=int, default=0, help='Start index in database')
p.add_argument('-i1', type=int, default=10000, help='Stop index in database')
p.add_argument('-bc', type=int, default=0, help='Boundary condition (neumann_zero=0, dirichlet_zero=1)')
args = p.parse_args()

for i, f in enumerate(args.filename):
    y = np.array([])
    q = np.array([])
    try:
        log = np.loadtxt(f"{f}_log.txt", skiprows=1)
    except IOError:
        break
    if (args.bc == 0):
        log_y = log[:,9]
    else:
        log_y = log[:,22]
    j = args.i0
    while (j >= 0):
        if (j > args.i1):
            break
        try:
            cross = np.loadtxt(f"{f}_cross_{str(j).zfill(6)}.txt", skiprows=1)
        except IOError:
            break
        cross = truncate_axis(cross, zmax=args.zmax, col=0)
        y = np.append(y, log_y[j])
        q = np.append(q, get_q_integral(cross))
        print(f"Processing {f}_cross_{str(j).zfill(6)}.txt")
        j += 1
    np.savetxt(f"y_q_integral_{f}.txt", np.column_stack([y, q]), header='y  q_integral', comments='', delimiter='  ', fmt='%.8e')    
    print(f"Saved y_q_integral_{f}.txt")
