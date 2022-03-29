#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 16:32:05 2021

@author: Baohong.Guo
"""

import numpy as np
import pandas as pd
import argparse
import abel
from scipy import interpolate

def get_ztop_index(Iz):
    for i in range(np.argmax(Iz), np.size(Iz)-2):
        if (Iz[i] < Iz[i+1]) & (Iz[i+1] < Iz[i+2]):
            index = i
            break
    else:
        index = np.size(Iz)-1
    return index

def get_FWHM_index(I):
    Iz = I.sum(axis=1)
    Ir = I.sum(axis=0)        
    z_index = np.amin(np.where(Iz >= (np.amax(Iz)+np.amin(Iz))/2))
    r_index = np.amax(np.where(Ir >= (np.amax(Ir)+np.amin(Ir))/2))
    return z_index, r_index

# =============================================================================
# def get_FWHM_index(I):
#     Ir = I.sum(axis=0)     
#     Ir0 = (np.amax(Ir)+np.amin(Ir))/2
#     Ir = np.amax(Ir)+np.amin(Ir)-Ir
#     r = np.arange(np.size(Ir))
#     i0 = np.amin(np.where(Ir >= (0.01*np.amax(Ir)+0.99*np.amin(Ir))))
#     i1 = np.amax(np.where(Ir <= (0.99*np.amax(Ir)+0.01*np.amin(Ir))))
#     r = r[i0:i1+1] 
#     Ir = Ir[i0:i1+1]
#     Iz = I.sum(axis=1)
#     Iz = Iz[0:np.argmax(Iz)+1]
#     Iz0 = (np.amax(Iz)+np.amin(Iz))/2
#     z = np.arange(np.size(Iz))
#     i0 = np.amin(np.where(Iz >= (0.001*np.amax(Iz)+0.999*np.amin(Iz))))
#     i1 = np.amax(np.where(Iz <= (0.999*np.amax(Iz)+0.001*np.amin(Iz))))
#     z = z[i0:i1+1]
#     Iz = Iz[i0:i1+1]
#     f_Ir = interpolate.UnivariateSpline(Ir, r, s=0)  
#     f_Iz = interpolate.UnivariateSpline(Iz, z, s=0)
#     r_index = f_Ir(Ir0)
#     z_index = f_Iz(Iz0)    
#     return z_index, r_index
# =============================================================================
    
p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("filename", type=str, nargs='+', help="Input filename(s)")
p.add_argument('-rmax', type=float, default=0.015, help='max(r) for plane')
p.add_argument('-zmax', type=float, default=0.125, help='max(z) for plane')
p.add_argument('-i0', type=int, default=0, help='Start index in database')
p.add_argument('-i1', type=int, default=10000, help='Stop index in database')
p.add_argument('-bc', type=int, default=0, help='Boundary condition (neumann_zero=0, dirichlet_zero=1)')
args = p.parse_args()

for i, f in enumerate(args.filename):
    t = np.array([])
    y_Emax = np.array([]) 
    R_Er = np.array([])
    y_3D = np.array([])
    R_3D = np.array([])
    y_2D = np.array([])  
    R_2D = np.array([])
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
            #I = np.loadtxt(f"{f}_plane_{str(j).zfill(6)}.vtk", dtype='str', skiprows=10)
            data = pd.read_csv(f"{f}_plane_{str(j).zfill(6)}.vtk", dtype='str', delim_whitespace=True, header=None, skiprows=10)
        except IOError:
            break
        data = data.dropna()
        data = np.array(data)
        I = data[0:int(np.shape(data)[0]/4)]
        for m in range(len(I)):
            try:
                I[m] = I[m].astype(float)
            except ValueError:
                for n in range(len(I[m])):
                    try:
                        I[m,n] = float(I[m,n])
                    except ValueError:
                        I[m,n] = 0
        dr = args.rmax/np.shape(I)[1]
        dz = args.zmax/np.shape(I)[0]
        Iz = I.sum(axis=1)
        I = I[:get_ztop_index(Iz)+1]
        t = np.append(t, log[j,1])
        y_Emax = np.append(y_Emax, log_y[j])
        R_Er = np.append(R_Er, log[j,14])
        if (j == 0):
            y_3D = np.append(y_3D, log_y[j])
            R_3D = np.append(R_3D, 0)
            y_2D = np.append(y_2D, log_y[j])
            R_2D = np.append(R_2D, 0)
        else:
            y_2D = np.append(y_2D, (get_FWHM_index(I)[0]+1)*dz)
            R_2D = np.append(R_2D, (get_FWHM_index(I)[1]+1)*dr)
            I_abel = abel.hansenlaw.hansenlaw_transform(I, direction='forward')
            y_3D = np.append(y_3D, (get_FWHM_index(I_abel)[0]+1)*dz)
            R_3D = np.append(R_3D, (get_FWHM_index(I_abel)[1]+1)*dr)
        print(f"Processing {f}_plane_{str(j).zfill(6)}.vtk ztop_index={get_ztop_index(Iz)} ztop={get_ztop_index(Iz)*dz}")
        j += 1
    np.savetxt(f"y_radius_N2_C3_{f}.txt", np.column_stack([t, y_Emax, R_Er, y_3D, R_3D, y_2D, R_2D]), 
               header='time  y_Emax  R_Er  y_3D  R_3D  y_2D  R_2D', comments='', delimiter='  ', fmt='%.8e')    
    print(f"Saved y_radius_N2_C3_{f}.txt")
