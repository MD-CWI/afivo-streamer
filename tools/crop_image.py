#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:33:37 2021

@author: Baohong.Guo
"""

import os   
import cv2 
 
# Define the function of CropImage
def CropImage(originpath, destpath, x0, x1, y0 ,y1):
    image = cv2.imread(originpath)  # Read the image
    name = originpath.split('\\')[-1]
    print(f'The shape of {name}:', image.shape)  # Print the shape of image
    CropImg = image[x0:x1, y0:y1]  # Crop the image to new shape
    cv2.imwrite(destpath, CropImg)  # Write the cropped image
           
if __name__ == '__main__':
    
    inputpath = r'C:\Users\Baohong.Guo\Desktop\1'  # The path of origin images
    outputpath = r'C:\Users\Baohong.Guo\Desktop\2'  # The path of cropped images
    
    # List the names of the entries in the given directory
    inputfile = os.listdir(inputpath)  
    for filename in inputfile:
        originpath = os.path.join(inputpath, filename)
        destpath = os.path.join(outputpath, filename)
        if os.path.isfile(originpath):
            CropImage(originpath, destpath, 0, 1024, 0, 1024)
        
