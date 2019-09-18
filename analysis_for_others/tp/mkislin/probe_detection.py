#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:54:58 2019

@author: wanglab
"""

import matplotlib as mpl, tifffile, numpy as np, cv2, matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter as gfilt
from scipy.ndimage import label
from tools.utils.io import load_kwargs
from tools.analysis.analyze_injection import orientation_crop_check, optimize_inj_detect

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

def find_site(im, thresh=3, filter_kernel=(3,3,3), num_sites_to_keep=1):

    if type(im) == str: im = tifffile.imread(im)

    filtered = gfilt(im, filter_kernel)
    thresholded = filtered > filtered.mean() + thresh*filtered.std() 
    labelled,nlab = label(thresholded)

    if nlab == 0:
        raise Exception('Site not detected, try a lower threshold?')
    elif nlab == 1:
        return labelled.astype(bool)
    elif num_sites_to_keep == 1:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        return labelled == np.argmax(sizes)+1
    else:
        sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
        vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
        return np.in1d(labelled, vals).reshape(labelled.shape)

src = "/home/wanglab/mounts/wang/oostland/lightsheet/m26/elastix/marlies_m26_1d3x_488_555_008na_1hfds_z5um_150msec_resized_ch01/result.tif"
src = orientation_crop_check(src, axes = ("2","0","1"), crop = "[:,423:,:]")
    
#optimize detection parameters for inj det
optimize_inj_detect(src, threshold=4, filter_kernel = (30,30,30))
#%%
brains = ["/jukebox/wang/oostland/lightsheet/m26","/jukebox/wang/oostland/lightsheet/m28"]

for brain in brains:
    kwargs = load_kwargs(brain)
    cellvol = [xx for xx in kwargs["volumes"] if xx.ch_type == "injch" or xx.ch_type == "cellch"][0]
    vol = tifffile.imread(cellvol.ch_to_reg_to_atlas)
    
    cb = vol[:, 423:, :]
