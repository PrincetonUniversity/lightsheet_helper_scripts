#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:04:02 2020

@author: wanglab
"""

import os, numpy as np, tifffile as tif, SimpleITK as sitk, cv2, multiprocessing as mp
from scipy.ndimage import zoom

# def resize_helper(img, dst, resizef):
#     print(os.path.basename(img))
#     im = sitk.GetArrayFromImage(sitk.ReadImage(img))
#     y,x = im.shape
#     yr = int(y/resizef); xr = int(x/resizef)
#     im = cv2.resize(im, (xr, yr), interpolation=cv2.INTER_LINEAR)
#     tif.imsave(os.path.join(dst, os.path.basename(img)), 
#                     im.astype("uint16"), compress=1)
    
# def load_memmap_arr(pth, mode="r", dtype = "uint16", shape = False):
#     """Function to load memmaped array.

#     Inputs
#     -----------
#     pth: path to array
#     mode: (defaults to r)
#     +------+-------------------------------------------------------------+
#     | "r"  | Open existing file for reading only.                        |
#     +------+-------------------------------------------------------------+
#     | "r+" | Open existing file for reading and writing.                 |
#     +------+-------------------------------------------------------------+
#     | "w+" | Create or overwrite existing file for reading and writing.  |
#     +------+-------------------------------------------------------------+
#     | "c"  | Copy-on-write: assignments affect data in memory, but       |
#     |      | changes are not saved to disk.  The file on disk is         |
#     |      | read-only.                                                  |
#     dtype: digit type
#     shape: (tuple) shape when initializing the memory map array

#     Returns
#     -----------
#     arr
#     """
#     if shape:
#         assert mode =="w+", "Do not pass a shape input into this function unless initializing a new array"
#         arr = np.lib.format.open_memmap(pth, dtype = dtype, mode = mode, shape = shape)
#     else:
#         arr = np.lib.format.open_memmap(pth, dtype = dtype, mode = mode)
#     return arr

# pth = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_488_Em_0/stitched/RES(7574x5773x3535)/102090/102090_120640"
# #path to store downsized images
# dst = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_488_Em_0/downsized"
# if not os.path.exists(dst): os.mkdir(dst)
# imgs = [os.path.join(pth, xx) for xx in os.listdir(pth)]
# z = len(imgs)
# resizef = 5 #factor to downsize imgs by
# iterlst = [(img, dst, resizef) for img in imgs]
# p = mp.Pool(12)
# p.starmap(resize_helper, iterlst)
# p.terminate()

#now downsample to 140% of allen atlas
dst = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_488_Em_0/downsized"
imgs = [os.path.join(dst, xx) for xx in os.listdir(dst)]; imgs.sort()
z = len(imgs)
y,x = sitk.GetArrayFromImage(sitk.ReadImage(imgs[0])).shape
arr = np.zeros((z,y,x))
atlpth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/average_template_25_sagittal_forDVscans.tif"
atl = sitk.GetArrayFromImage(sitk.ReadImage(atlpth))
atlz,atly,atlx = atl.shape #get shape, sagittal
#read all the downsized images
for i,img in enumerate(imgs):
    if i%10==0: print(i)
    arr[i,:,:] = sitk.GetArrayFromImage(sitk.ReadImage(img)) #horizontal
#switch to sagittal
arr = np.swapaxes(arr,2,0)
print("\n**********downsizing....heavy!**********\n")

arrd = zoom(arr, ((atlz*1.4/z),(atly*1.4/y),(atlx*1.4/x)), order=1)
tif.imsave(os.path.join(os.path.dirname(dst), "downsized_for_atlas.tif"), arrd.astype("uint16"))
