#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 12:51:13 2019

@author: wanglab
"""
from skimage.external import tifffile
import numpy as np

pth = "/home/wanglab/Documents/electrode.tif"

img = tifffile.imread(pth).astype("int16")

zr,yr,xr = np.where(img < 0)

for z,y,x in zip(zr,yr,xr):
    img[z,y,x] = 65000

            
tifffile.imsave("/home/wanglab/Documents/electrode_corrected.tif", img.astype("uint16"))