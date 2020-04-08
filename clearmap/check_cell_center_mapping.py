#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:09:29 2018

@author: wanglab
"""

import numpy as np, os, time, cv2, pandas as pd
from skimage.external import tifffile
import matplotlib.pyplot as plt

def resize_merged_stack(pth, dst, dtype = "uint16", resizef = 3):
    """
    resize function for large image stacks using cv2
    inputs:
        pth = 4d stack, memmap array or numpy array
        dst = path of tif file to save
        dtype = default uint16
        resizef = default 6
    """

    #read file
    if pth[-4:] == ".tif": img = tifffile.imread(pth)
    elif pth[-4:] == ".npy": img = np.lib.format.open_memmap(pth, dtype = dtype, mode = "r")
    else: img = pth #if array was input

    z,y,x,ch = img.shape
    resz_img = np.zeros((z, int(y/resizef), int(x/resizef), ch))

    for i in range(z):
        for j in range(ch):
            #make the factors - have to resize both image and cell center array
            xr = int(img[i, :, :, j].shape[1] / resizef); yr =  int(img[i, :, :, j].shape[0] / resizef)
            im = cv2.resize(img[i, :, :, j], (xr, yr), interpolation=cv2.INTER_LINEAR)
            resz_img[i, :, :, j] = im.astype(dtype)

    tifffile.imsave(dst, resz_img.astype(dtype))

    return dst

def check_cell_center_to_fullsizedata(brain, zstart, zstop, dst, resizef):
    """
    maps cnn cell center coordinates to full size cell channel images
    inputs:
        brain = path to lightsheet processed directory
        zstart = beginning of zslice
        zstop = end of zslice
        dst = path of tif stack to save
    NOTE: 20+ PLANES CAN OVERLOAD MEMORY
    """
    start = time.time()

    #doing things without loading parameter dict
    fzfld = os.path.join(brain, "full_sizedatafld")

    #exception if only 1 channel is imaged
    cellch = os.path.join(fzfld, [xx for xx in os.listdir(fzfld) if "647" in xx][0])

    #not the greatest way to do things, but works
    src = [os.path.join(cellch, xx) for xx in os.listdir(cellch) if xx[-3:] == "tif" and int(xx[-7:-4]) in range(zstart, zstop)]; src.sort()

    raw = np.zeros((len(src), tifffile.imread(src[0]).shape[0], tifffile.imread(src[0]).shape[1]))

    for i in range(len(src)):
        raw[i, :, :] = tifffile.imread(src[i])

    pth = os.path.join(brain, "clearmap_cluster_output/cells.npy")
    cells = np.load(pth) #this is in x,y,z!!!!!!!!!!!!!!!

    cells = cells[(cells[:, 2] >= zstart) & (cells[:, 2] <= zstop-1)] #-1 to account for range

    cell_centers = np.zeros(raw.shape)

    for i, r in enumerate(cells):
        cell_centers[r[2]-zstart, r[1]-1:r[1]+1, r[0]-1:r[0]+1] = 50000

    rbg = np.stack([raw.astype("uint16"), cell_centers.astype("uint16"), np.zeros_like(raw)], -1)

    resize_merged_stack(rbg, os.path.join(dst, "{}_raw_cell_centers_resizedfactor{}_z{}-{}.tif".format(os.path.basename(brain),
                                          resizef, zstart, zstop)), "uint16", resizef)

    print("took %0.1f seconds to make merged maps for %s" % ((time.time()-start), brain))

def check_cell_center_to_resampled(brain, zstart, zstop, dst):
    """
    maps cnn cell center coordinates to resampled stack
    inputs:
        brain = path to lightsheet processed directory
        zstart = beginning of zslice
        zstop = end of zslice
        dst = path of tif stack to save
    NOTE: 20+ PLANES CAN OVERLOAD MEMORY
    """
    start = time.time()

    #doing things without loading parameter dict, could become a problem
    tifs = [xx for xx in os.listdir(brain) if xx[-4:] == ".tif"]; tifs.sort()
    raw = tifffile.imread(tifs[len(tifs)-1])

    pth = os.path.join(brain, "clearmap_cluster_output/cells.npy")
    cells = np.load(pth) #this is in x,y,z!!!!!!!!!!!!!!!

    cells = cells[(cells[:, 2] >= zstart) & (cells[:, 2] <= zstop-1)] #-1 to account for range

    cell_centers = np.zeros(raw.shape)

    for i, r in enumerate(cells):
        cell_centers[r[2]-zstart, r[1]-1:r[1]+1, r[0]-1:r[0]+1] = 50000

    rbg = np.stack([raw.astype("uint16"), cell_centers.astype("uint16"), np.zeros_like(raw)], -1)

    resize_merged_stack(rbg, os.path.join(dst, "{}_raw_cell_centers_resized_z{}-{}.tif".format(os.path.basename(brain),
                                          zstart, zstop)), "uint16", 6)

    print("took %0.1f seconds to make merged maps for %s" % ((time.time()-start), brain))


if __name__ == "__main__":
    src = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/processed/"
    dst = "/jukebox/wang/Jess/lightsheet_output/201904_ymaze_cfos/qc"
    if not os.path.exists(dst): os.mkdir(dst)
    zstart = 300; zstop = 320 #z-plane #'s you want to visualize at a time, 250-400 is probably ideal for thalamus
    resizef = 1 #factor by which to downsize raw and cell center overlay, keep 1 if you do not want to downsize
    ids = ["an20"]

    for i in ids:
        brain = os.path.join(src, i)
        check_cell_center_to_fullsizedata(brain, zstart, zstop, dst, resizef)
