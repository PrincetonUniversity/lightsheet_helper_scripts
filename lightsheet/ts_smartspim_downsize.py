#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:04:02 2020

@author: wanglab
"""

import os, numpy as np, tifffile as tif, SimpleITK as sitk, cv2, multiprocessing as mp

def resize_helper(img, dst, resizef):
    print(os.path.basename(img))
    im = sitk.GetArrayFromImage(sitk.ReadImage(img))
    y,x = im.shape
    yr = int(y/resizef); xr = int(x/resizef)
    im = cv2.resize(im, (xr, yr), interpolation=cv2.INTER_LINEAR)
    tif.imsave(os.path.join(dst, os.path.basename(img)), 
                    im.astype("uint16"), compress=1)
    
def load_memmap_arr(pth, mode="r", dtype = "uint16", shape = False):
    """Function to load memmaped array.

    Inputs
    -----------
    pth: path to array
    mode: (defaults to r)
    +------+-------------------------------------------------------------+
    | "r"  | Open existing file for reading only.                        |
    +------+-------------------------------------------------------------+
    | "r+" | Open existing file for reading and writing.                 |
    +------+-------------------------------------------------------------+
    | "w+" | Create or overwrite existing file for reading and writing.  |
    +------+-------------------------------------------------------------+
    | "c"  | Copy-on-write: assignments affect data in memory, but       |
    |      | changes are not saved to disk.  The file on disk is         |
    |      | read-only.                                                  |
    dtype: digit type
    shape: (tuple) shape when initializing the memory map array

    Returns
    -----------
    arr
    """
    if shape:
        assert mode =="w+", "Do not pass a shape input into this function unless initializing a new array"
        arr = np.lib.format.open_memmap(pth, dtype = dtype, mode = mode, shape = shape)
    else:
        arr = np.lib.format.open_memmap(pth, dtype = dtype, mode = mode)
    return arr

pth = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_488_Em_0/stitched/RES(7574x5773x3535)/102090/102090_120640"
#path to store downsized images
dst = "/jukebox/LightSheetTransfer/kelly/2020_07_15/20200715_12_14_06_f37080_mouse2_20171015/Ex_488_Em_0/downsized"
if not os.path.exists(dst): os.mkdir(dst)
imgs = [os.path.join(pth, xx) for xx in os.listdir(pth)]
z = len(imgs)
resizef = 5 #factor to downsize imgs by
iterlst = [(img, dst, resizef) for img in imgs]
p = mp.Pool(12)
p.starmap(resize_helper, iterlst)
p.terminate()
