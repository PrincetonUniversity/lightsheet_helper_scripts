#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:09:29 2018

@author: wanglab
"""

import numpy as np, os, time, cv2
from skimage.external import tifffile
import matplotlib.pyplot as plt
import pandas as pd

def resize_merged_stack(pth, dst, dtype = "uint16", resizef = 6):
    """ 
    resize function for large image stacks using cv2
    inputs:
        pth = 3d stack or memmap array or numpy array
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
            #make the factors
            xr = int(img[i, :, :, j].shape[1] / resizef); yr =  int(img[i, :, :, j].shape[0] / resizef)
            im = cv2.resize(img[i, :, :, j], (xr, yr), interpolation=cv2.INTER_LINEAR)
            resz_img[i, :, :, j] = im.astype(dtype)
    
    tifffile.imsave(dst, resz_img.astype(dtype))
    
    return dst

def check_cell_center_to_fullsizedata(brain, zstart, zstop, dst):
    """ 
    maps cnn cell center coordinates to full size cell channel images
    inputs:
        brain = path to lightsheet processed directory
        zstart = beginning of zslice
        zstop = end of zslice
        dst = path of tif stack to save
    NOTE: ONLY WORKS FOR 5 Z PLANES AT A TIME (#FIXME)
    """
    start = time.time()
    
    data = os.listdir(os.path.join(brain, "full_sizedatafld")); data.sort()
    cellch = os.path.join(brain, "full_sizedatafld/"+data[2])
    
    #really messy way to do things, also only limited to 5 planes    
    src = [os.path.join(cellch, xx) for xx in os.listdir(cellch) if xx[-9:-4] == "Z0"+str(zstart) or xx[-9:-4] == "Z0"+str(zstart+1) or
           xx[-9:-4] == "Z0"+str(zstart+2) or xx[-9:-4] == "Z0"+str(zstart+3) or xx[-9:-4] == "Z0"+str(zstop)]
    
    raw = np.zeros((len(src), tifffile.imread(src[0]).shape[0], tifffile.imread(src[0]).shape[1]))
    
    for i in range(len(src)):
        raw[i, :, :] = tifffile.imread(src[i])
        
    pth = os.path.join(brain, "3dunet_output/pooled_cell_measures/"+os.path.basename(brain)+"_cell_measures.csv")
    cells = pd.read_csv(pth)
    
    plt.imshow(raw[0,:,:])
    cells = cells[(cells["z"] >= zstart) & (cells["z"] <= zstop)]
    
    cell_centers = np.zeros(raw.shape)
    
    for i, r in cells.iterrows():
        cell_centers[r["z"]-zstart, r["y"]-5:r["y"]+5, r["x"]-5:r["x"]+5] = 1
        
    rbg = np.stack([raw.astype("uint16"), cell_centers.astype("uint16"), np.zeros_like(raw)], -1)

    resize_merged_stack(rbg, os.path.join(dst, "{}_raw_cell_centers_resized_z{}-{}.tif".format(os.path.basename(brain), 
                                          zstart, zstop)), "uint16", 6)
    
    print("took {} seconds to make merged maps".format(time.time() - start))
    
if __name__ == "__main__":

    brain = "/jukebox/wang/pisano/tracing_output/antero_4x/20180612_jg77"
    zstart = 400; zstop = 404
    dst = "/home/wanglab/Desktop"
    
    check_cell_center_to_fullsizedata(brain, zstart, zstop, dst)