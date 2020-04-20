#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 13:29:33 2019

@author: wanglab
"""

import os, tifffile, matplotlib.pyplot as plt, numpy as np, cv2
import matplotlib.colors
from skimage.morphology import ball

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "lime"]) #lime color makes cells pop

pth = "/jukebox/wang/zahra/h129_qc/thal_transformed_points"

brains_to_inspect = ["20180410_jg51_bl6_lob6b_04", "20180417_jg59_bl6_cri_03"]

brain_pths = [os.path.join(os.path.join(pth, xx), "transformed_volume/merged.tif") for xx in brains_to_inspect]

i=1
name = os.path.basename(os.path.dirname(os.path.dirname(brain_pths[i])))
brain = tifffile.imread(brain_pths[i])
img = brain[200:204,:,:,1]
cell = brain[200:204,:,:,0]

#apply x y dilation
r = 2
selem = ball(r)[int(r/2)]
cell = cell.astype("uint8")
cell = np.asarray([cv2.dilate(cell[i], selem, iterations = 1) for i in range(cell.shape[0])])

merged = np.stack([np.max(img, axis=0), np.max(cell, axis=0), np.zeros_like(np.max(cell, axis=0))], -1)

tifffile.imsave("/home/wanglab/Desktop/{}.tif".format(name), merged)
