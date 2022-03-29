#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 18:45:07 2021

@author: Baohong.Guo
"""

import os
import cv2
import numpy as np
from scipy import integrate
from scipy import interpolate
from scipy.signal import savgol_filter


# The log array should be in descending order, otherwise reversing it using log[::-1] or np.flipud(log)
# Get the max index of descend array based on Emax
def get_descend_Emax_max_index(Emax, Emax_min=3.5e6):  
    index = np.size(Emax)-1
    for i in range(np.size(Emax)):
        if (Emax[i] < Emax_min):
            index = i-1
            break
    return index   # Return the max index

# Get the max index of descend array based on y
def get_descend_max_index(y, ymin=0.03): 
    index1 = np.size(y)-1
    for i in range(np.size(y)):
        if (i < np.size(y)-1):
            if (y[i] < y[i+1]) & (y[i] > 5e-3):             
                index1 = i
                break
    i1 = np.amax(np.append(np.where(y >= ymin), 0))
    i2 = np.amin(np.append(np.where(y <= ymin), np.size(y)-1))
    if np.abs(y[i1]-ymin) <= np.abs(y[i2]-ymin):
        index2 = i1
    else:
        index2 = i2
    index = np.minimum(index1, index2)
    return index   # Return the max index

# Get the min index of descend array
def get_descend_min_index(y, ymax=0.1151):
    i1 = np.amax(np.append(np.where(y >= ymax), 0))
    i2 = np.amin(np.append(np.where(y <= ymax), np.size(y)-1))
    if np.abs(y[i1]-ymax) <= np.abs(y[i2]-ymax):
        index = i1
    else:
        index = i2
    return index   # Return the min index

# Truncate log array based on y
def truncate_log(log, ymax=1e3, ymin=0.03, col=9): 
    log = log.astype(float)
    y = log[:,int(col)]
    i1 = get_descend_min_index(y, ymax)  
    i2 = get_descend_max_index(y, ymin) 
    log = log[i1:i2+1]  
    return log  # Return truncated log

# Truncate log array based on Emax
def truncate_Emax_log(log, ymax=1e3, ymin=0.03, Emax_min=3.5e6, col=9, col_Emax=7): 
    log = log.astype(float)
    y = log[:,int(col)]
    Emax = log[:,int(col_Emax)]
    i1 = get_descend_min_index(y, ymax)  
    i2 = get_descend_max_index(y, ymin) 
    i3 = get_descend_Emax_max_index(Emax, Emax_min)
    log = log[i1:np.minimum(i2,i3)+1]  
    return log  # Return truncated log

# Get the time information of max index
def get_time(log, ymin=0.03, Emax_min=3.5e6, col=9, col_Emax=7): 
    log = log.astype(float)
    y = log[:,int(col)]
    Emax = log[:,int(col_Emax)]
    i1 = get_descend_max_index(y, ymin)
    i2 = get_descend_Emax_max_index(Emax, Emax_min)
    i = np.minimum(i1,i2)
    it = log[:,0].astype(int)
    time = log[:,1]
    v_savgol = np.abs(savgol_filter(y, 11, 2, deriv=1, delta=time[1]-time[0]))
    return it[i], time[i], y[i], Emax[i], v_savgol[i]  # Return the time information

# The z array should be in ascending order, otherwise reversing it using z[::-1] or np.flipud(z)
# Get the nearest index of ascend array
def get_ascend_index(z, z_value=0.03):
    z = z.astype(float)      
    i1 = np.amax(np.append(np.where(z <= z_value), 0))
    i2 = np.amin(np.append(np.where(z >= z_value), np.size(z)-1))
    if np.abs(z[i1]-z_value) <= np.abs(z[i2]-z_value):
        index = i1
    else:
        index = i2
    return index

# Truncate axis array
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
    return array[i1:i2] # Return truncated array

# Calculate the radius from the electrc field ahead of streamer
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

# Integrate line charge density along streamer axis for a couple of timestep
def integrate_λ(crossfilepath, logfile, zmax=0.1151, zmin=0):
    y = np.array([])
    q = np.array([])
    log = np.loadtxt(logfile, skiprows=1)
    crossfile = os.listdir(crossfilepath)
    for i, f in enumerate(crossfile):
        filename = os.path.join(crossfilepath, f)
        cross = np.loadtxt(filename, skiprows=1)
        cross = truncate_axis(cross, zmax, zmin, col=0)
        z = cross[:,0]
        λ = cross[:,2]*1.6e-19
        q = np.append(q, integrate.trapz(λ, z))
        y = np.append(y, log[i,9])
    q_integrate = np.column_stack([y, q])
    return q_integrate

# Integrate line charge density at a moment
def get_q_integral(cross):
    z = cross[:,0]
    λ = cross[:,2]*1.6e-19
    q = integrate.trapz(λ, z)
    return q

# Get the index of FWHM
def get_FWHM_index(I):
    Iz = I.sum(axis=1)
    Ir = I.sum(axis=0)        
    z_index = np.amin(np.where(Iz >= (np.amax(Iz)+np.amin(Iz))/2))
    r_index = np.amax(np.where(Ir >= (np.amax(Ir)+np.amin(Ir))/2))
    return z_index, r_index

# Crop image
def crop_image(originpath, destpath, x0=0, x1=100, y0=0 ,y1=100): 
    image = cv2.imread(originpath)  # Read the image 
    imagename = originpath.split('\\')[-1]
    print(f'The shape of {imagename}:', image.shape)  # Print the shape of image
    cropimg = image[x0:x1, y0:y1]  # Crop the image into new shape
    cv2.imwrite(destpath, cropimg)  # Write the cropped image

# Generate electric field profile versus y
def generate_y_field(profile='exponential', initial_field=1.15e6, goal_field=0.2932e6, decay_distance=1.1e-3, z_constant=1.2e-3, z_fall=46.6e-3):
    npoint = 5000
    z = np.linspace(z_constant + z_fall, 0, npoint)
    l = np.maximum(z_fall-z, 0)
    field = np.zeros(npoint)
    if (profile == 'exponential'):    
        field = goal_field + (initial_field - goal_field) * np.exp(-l/decay_distance)
    elif (profile == 'linear'):
        field = initial_field + decay_distance * l
        j = npoint
        for i in range(npoint):
            if ((initial_field >= goal_field) & (field[i] <= goal_field)):
                j=i; break
            if ((initial_field <= goal_field) & (field[i] >= goal_field)):
                j=i; break
        for i in range(j, npoint):
            field[i] = goal_field          
    elif (profile == 'step'):
        for i in range(npoint):
            if (l[i] == 0):
                field[i] = initial_field
            else:
                field[i] = goal_field
    else:  # Return constant field for other profiles
        for i in range(npoint):
            field[i] = initial_field
    return np.column_stack([z, field])  # Return array containing z and field

# Generate electric field profile versus time
def generate_t_field(tStart=5.0e-9, tRise=5.0e-9, tPulse=10.0e-9, tFall=5.0e-9, tPost=10.0e-9, \
                     field_start=1.0e6, field_pulse=2.0e6, field_post=0.0e6):
    nPoints = 5000
    tTotal = tStart + tRise + tPulse + tFall + tPost
    time = np.linspace(0.0, tTotal, nPoints)
    field = np.zeros(nPoints)
    for i in range(nPoints):
        if (time[i] <= tStart):
          field[i] = field_start
        elif (time[i] > tStart and time[i] <= (tStart + tRise)):
          field[i] = field_start + (field_pulse - field_start) / tRise * (time[i] - tStart)
        elif(time[i] > (tStart + tRise) and time[i] <= (tStart + tRise + tPulse)):
          field[i] = field_pulse
        elif(time[i] > (tStart + tRise + tPulse) and time[i] <= (tStart + tRise + tPulse + tFall)):
          field[i] = field_pulse + (field_post - field_pulse) / tFall * (time[i] - (tStart + tRise + tPulse))
        else:
          field[i] = field_post
    return np.column_stack([time, field])