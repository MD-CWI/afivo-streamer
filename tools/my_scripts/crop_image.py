#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:33:37 2021

@author: Baohong.Guo
"""

import os   
from my_module import crop_image
        
if __name__ == '__main__':
    
    inputpath = r'visit_save'  # The path of origin images
    outputpath = r'crop'  # The path of cropped images
    
    # List the names of the entries in the given directory
    inputfile = os.listdir(inputpath)  
    for filename in inputfile:
        originpath = os.path.join(inputpath, filename)
        destpath = os.path.join(outputpath, filename)
        if os.path.isfile(originpath):
            crop_image(originpath, destpath, 46, 806, 195, 320)