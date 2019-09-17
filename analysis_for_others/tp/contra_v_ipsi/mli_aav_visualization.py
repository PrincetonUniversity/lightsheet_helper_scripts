#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:03:54 2019

@author: wanglab
"""

import os, tifffile, numpy as np, matplotlib.pyplot as plt
from tools.utils.io import load_kwargs

src = "/jukebox/wang/Jess/lightsheet_output/201906_development_cno/processed"
dst = "/home/wanglab/Desktop"

animals = ["an01", "an02", "an03", "an04", "an05", "an06", "an07", "an09", "an10", "an12", "an13", "an14", "an15", "an16", "an17"]

brains = [os.path.join(src, xx) for xx in animals]

#coronal
maxip_scale = 4 #aka 60 um sections
#thalamus
maxip_start = 350
maxip_stop = 450

cb_start = 600
cb_stop = 700

for brain in brains:
    
    kwargs = load_kwargs(brain)
    cellvol = [xx for xx in kwargs["volumes"] if xx.ch_type == "injch" or xx.ch_type == "cellch"][0]
    vol = tifffile.imread(cellvol.resampled_for_elastix_vol)
    vol_cor = np.rot90(np.transpose(vol, [1, 0, 2]), axes = (2,1))
    
    #init figure
    f = np.sqrt((maxip_stop - maxip_start)/maxip_scale).astype(int) #has to be a square 
    ncols, nrows = f, f
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = (15,10), sharex = True, gridspec_kw = {"wspace":0, "hspace":0})
    slcs = np.arange(maxip_start, maxip_stop, maxip_scale)
    k = 0 #init slice range
    for i in range(ncols):
        for j in range(nrows):
            axes[i,j].pcolormesh(np.max(vol_cor[slcs[k]:slcs[k]+maxip_scale]*25, axis=0), cmap="Greys")
            axes[i,j].invert_yaxis()
            axes[i,j].axis("off")
            k += 1        
    
    #done with the page
    plt.savefig(os.path.join(dst, os.path.basename(brain)+"_aav_thalz%d-%d_zstep%d.png" % (maxip_start, 
                                      maxip_stop, maxip_scale)), bbox_inches = "tight")             
    
    plt.close()
    
    #init figure
    f = np.sqrt((cb_stop - cb_start)/maxip_scale).astype(int) #has to be a square 
    ncols, nrows = f, f
    fig, axes = plt.subplots(ncols = ncols, nrows = nrows, figsize = (15,10), sharex = True, gridspec_kw = {"wspace":0, "hspace":0})
    slcs = np.arange(cb_start, cb_stop, maxip_scale)
    k = 0 #init slice range
    for i in range(ncols):
        for j in range(nrows):
            axes[i,j].pcolormesh(np.max(vol_cor[slcs[k]:slcs[k]+maxip_scale]*15, axis=0), cmap="Greys")
            axes[i,j].invert_yaxis()
            axes[i,j].axis("off")
            k += 1        
    
    #done with the page
    plt.savefig(os.path.join(dst, os.path.basename(brain)+"_aav_cbz%d-%d_zstep%d.png" % (cb_start, cb_stop, 
                             maxip_scale)), bbox_inches = "tight")             
    plt.close()
    