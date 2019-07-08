#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:47:05 2019

@author: wanglab
"""

import os, numpy as np, tifffile, SimpleITK as sitk, matplotlib.pyplot as plt
%matplotlib inline

#path to raw data
pth = "/jukebox/wang/seagravesk/lightsheet/cfos_raw_images/cfos/171206_f37077_observer_20171011_790_015na_1hfsds_z5um_1000msec_12-27-06"

#find only raw images
imgs = [os.path.join(pth, xx) for xx in os.listdir(pth) if "UltraII_raw_DataStack" in xx 
        or "UltraII_raw_RawDataStack" in xx]; imgs.sort()

y,x = tifffile.imread(imgs[2]).shape
print(y,x)
#make memmap arr of these
#only left ls, one channel, one tile
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/bigstitcher/f37077_observer_20171011_dorsal_up_790.npy", 
                                dtype = "uint16", mode = "w+", shape = (len(imgs), y, x))

for i, img in enumerate(imgs):
    #read 1 plane
    im = sitk.GetArrayFromImage(sitk.ReadImage((img)))
    arr[i] = im
    print(i)
    if i%50 == 0: arr.flush()
    
#%%
import h5py 
    
#read mmap array and make into raw file?
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/bigstitcher/f37077_observer_20171011_dorsal_up_790.npy", 
                                dtype = "uint16", mode = "r")

#h = h5py.File("/jukebox/scratch/zmd/bigstitcher/f37077_observer_20171011_dorsal_up_790.hdf5", 'w')   
#dset = h.create_dataset("data", data = arr)

tifffile.imsave("/jukebox/scratch/zmd/bigstitcher/f37077_observer_20171011_dorsal_up_790_z500_600.tif", arr[500:600])

#%%
#read mmap array and make into raw file?
arr = np.lib.format.open_memmap("/jukebox/scratch/zmd/bigstitcher/f37077_observer_20171011_ventral_up_790.npy", 
                                dtype = "uint16", mode = "r")

#h = h5py.File("/jukebox/scratch/zmd/bigstitcher/f37077_observer_20171011_ventral_up_790.hdf5", 'w')   
#dset = h.create_dataset("data", data = arr)
arr_flip = np.flipud(np.fliplr(arr))
tifffile.imsave("/jukebox/scratch/zmd/bigstitcher/f37077_observer_20171011_ventral_up_790.tif", arr_flip)
    