#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:47:05 2019

@author: wanglab
"""

import os, numpy as np, tifffile, SimpleITK as sitk, matplotlib.pyplot as plt
%matplotlib inline

#path to raw data
pth = "/jukebox/LightSheetData/mallarino/ricardo/20190702/raw_data/190702_20190313_edu_171_20190610_mallarino_1d3x_488_008na_z25um_1hfds_100ms_40povlp_sagittal_11-03-29"

#find only raw images for tile 0
tile0 = [os.path.join(pth, xx) for xx in os.listdir(pth) if "UltraII_raw_RawDataStack[00 x 00]_C01" in xx]; tile0.sort()
print("\n\ntile 0: {} images".format(len(tile0)))
y,x = tifffile.imread(tile0[2]).shape
print(y,x)
#make memmap arr of these
#only left ls, one channel, one tile
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/bigstitcher/glider/20190313_edu_171_488_tile00_c01.npy", 
                                dtype = "uint16", mode = "w+", shape = (len(tile0), y, x))

for i, img in enumerate(tile0):
    #read 1 plane
    im = sitk.GetArrayFromImage(sitk.ReadImage((img)))
    arr[i] = im
    print(i)
    if i%50 == 0: arr.flush()    
    
#find only raw images for tile 1
tile1 = [os.path.join(pth, xx) for xx in os.listdir(pth) if "UltraII_raw_RawDataStack[01 x 00]_C01" in xx]; tile1.sort()
print("\n\ntile 1: {} images".format(len(tile1)))
y,x = tifffile.imread(tile1[2]).shape
print(y,x)
#make memmap arr of these
#only left ls, one channel, one tile
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/bigstitcher/glider/20190313_edu_171_488_tile01_c01.npy", 
                               dtype = "uint16", mode = "w+", shape = (len(tile1), y, x))
for i, img in enumerate(tile1):
    #read 1 plane
    im = sitk.GetArrayFromImage(sitk.ReadImage((img)))
    arr[i] = im
    print(i)
    if i%50 == 0: arr.flush()

#%%

#read mmap array and make into raw file?
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/bigstitcher/glider/20190313_edu_171_488_tile00_c01.npy", 
                                dtype = "uint16", mode = "r")

tifffile.imsave("/jukebox/scratch/zmd/bigstitcher/glider/20190313_edu_171_488_tile00_c01.tif", arr)

#read mmap array and make into raw file?
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/bigstitcher/glider/20190313_edu_171_488_tile01_c01.npy", 
                                dtype = "uint16", mode = "r")

tifffile.imsave("/jukebox/scratch/zmd/bigstitcher/glider/20190313_edu_171_488_tile01_c01.tif", arr)
    