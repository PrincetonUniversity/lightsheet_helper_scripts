#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:47:05 2019

@author: wanglab
"""

import os, numpy as np, tifffile, SimpleITK as sitk, matplotlib.pyplot as plt
%matplotlib inline

#path to raw data
pth = "/jukebox/LightSheetData/brodyatlas/raw_data/190613_a238_201904_1d3x_488_008na_1hfds_z10um_50msec_20povlp_11-01-35"

#bad way to do this, use regex
#find only raw images for tile 0
tile0 = [os.path.join(pth, xx) for xx in os.listdir(pth) if "UltraII_raw_RawDataStack[00 x 00]_C01" in xx]; tile0.sort()
print("\n\ntile 0: {} images".format(len(tile0)))

#checks
img0 = tifffile.imread(tile0[300])
plt.imshow(img0)
y,x = img0.shape
print(y,x)

#make memmap arr of these
#only left ls, one channel, one tile
dst0 = "/jukebox/scratch/zmd/bigstitcher/rat/a238_201904_tile00_c01.npy"
arr = np.lib.format.open_memmap(dst0, dtype = "uint16", mode = "w+", shape = (len(tile0), y, x))

for i, img in enumerate(tile0):
    #read 1 plane
    im = sitk.GetArrayFromImage(sitk.ReadImage((img)))
    arr[i] = im
    print(i)
    if i%200 == 0: arr.flush()    
    
#find only raw images for tile 1
tile1 = [os.path.join(pth, xx) for xx in os.listdir(pth) if "UltraII_raw_RawDataStack[01 x 00]_C01" in xx]; tile1.sort()
print("\n\ntile 1: {} images".format(len(tile1)))

#checks
img1 = tifffile.imread(tile1[300])
plt.imshow(img1)
y,x = img1.shape
print(y,x)

#make memmap arr of these
#only left ls, one channel, one tile
dst1 = "/jukebox/scratch/zmd/bigstitcher/rat/a238_201904_tile01_c01.npy"
arr = np.lib.format.open_memmap(dst1, dtype = "uint16", mode = "w+", shape = (len(tile1), y, x))
for i, img in enumerate(tile1):
    #read 1 plane
    im = sitk.GetArrayFromImage(sitk.ReadImage((img)))
    arr[i] = im
    print(i)
    if i%200 == 0: arr.flush()

#%%
#read mmap array and make into raw file?
arr = np.lib.format.open_memmap(dst0, dtype = "uint16", mode = "r")

tifffile.imsave(dst0[:-3]+"tif", arr)

#read mmap array and make into raw file?
arr = np.lib.format.open_memmap(dst1, dtype = "uint16", mode = "r")

tifffile.imsave(dst1[:-3]+"tif", arr)

#%%

#playing around with metadata
test = "/home/wanglab/mounts/LightSheetData/brodyatlas/raw_data/190613_a238_201904_1d3x_488_008na_1hfds_z10um_50msec_20povlp_11-01-35/11-01-35_UltraII_raw_RawDataStack[00 x 00]_C00_xyz-Table Z0000_UltraII Filter0000.ome.tif"
with tifffile.TiffFile(test) as tif:
    metadata = tif.ome_metadata
#    metadata = tif[0].image_description
#metadata = json.loads(metadata.decode('utf-8'))
