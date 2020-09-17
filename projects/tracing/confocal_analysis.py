#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 16:09:32 2020

@author: wanglab
"""

import os, SimpleITK as sitk, numpy as np, tifffile, matplotlib.pyplot as plt, multiprocessing as mp

pth = "/home/wanglab/wang/zahra/confocal/slice6"

#this is without blending though
strngs = np.unique(np.array([xx[-11:] for xx in os.listdir(pth) if "tif" in xx]))

for strng in strngs:
    print(strng)
    imgs = [os.path.join(pth,xx) for xx in os.listdir(pth) if strng in xx]; imgs.sort()
    img = sitk.GetArrayFromImage(sitk.ReadImage(imgs[0]))
    #init empty stitched image
    arr = np.zeros((img.shape[0]*16, img.shape[0]*20)) #y,x order
    
    # counter for images
    i=0
    for col in range(16):
        print(col)
        for row in range(20):
            arr[int(col*512):int((col+1)*512),int(row*512):int((row+1)*512)] = sitk.GetArrayFromImage(sitk.ReadImage(imgs[i]))
            i+=1
    
    tifffile.imsave("/home/wanglab/Desktop/freeksimo_m1-16-14838-02_2017-01-15_slice6_lefthalf_merged_{}".format(strng),arr.astype("uint8"))

#%%
# read and resave all images?
def save(impth):
    im = sitk.GetArrayFromImage(sitk.ReadImage(impth))
    src = "/home/wanglab/Desktop/slice6"
    tifffile.imsave(os.path.join(src, os.path.basename(impth)), im.astype("uint8").T) #swapping x y!!
    
    return
pth = "/home/wanglab/wang/zahra/confocal/slice6"
imgs = [os.path.join(pth,xx) for xx in os.listdir(pth) if "tif" in xx]; imgs.sort()

p = mp.Pool(12)
p.map(save,imgs)