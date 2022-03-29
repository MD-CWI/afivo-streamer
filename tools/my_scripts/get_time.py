#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:32:05 2021

@author: Baohong.Guo
"""

import numpy as np
import os
from my_module import get_time

if __name__ == '__main__':
    
    filepath = r'../data_files/log_files'  # The path of files 
    
    # List the names of the entries in the given directory
    allfile = os.listdir(filepath)
    print('filename  it  time  y  Emax  v_savgol')
    for f in allfile:
        filename = os.path.join(filepath, f)
        log = np.loadtxt(filename, skiprows=1)
        print(f, get_time(log, ymin=0.1, Emax_min=3.5e6, col=9, col_Emax=7))