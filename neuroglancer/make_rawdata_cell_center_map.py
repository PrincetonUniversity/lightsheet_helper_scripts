#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:51:53 2020

@author: wanglab
"""

import numpy as np, os, pandas as pd
from skimage.external import tifffile
from scipy.ndimage import gaussian_filter

src = "/jukebox/wang/pisano/tracing_output/antero_4x"
dst = "/jukebox/scratch/zmd"
brains = ["20180410_jg51_bl6_lob6b_04"]

#for array job parallelization
print(os.environ["SLURM_ARRAY_TASK_ID"])
jobid = int(os.environ["SLURM_ARRAY_TASK_ID"])

brain = brains[jobid]
print("\n"+brain+"\n")
#join paths
brain_pth = os.path.join(src, brain)

#cell dimensions
cells = pd.read_csv(os.path.join(brain_pth, "3dunet_output/pooled_cell_measures/%s_cell_measures.csv" % brain))
fullszfld = os.path.join(brain_pth, "full_sizedatafld")
cellfld = [os.path.join(fullszfld, xx) for xx in os.listdir(fullszfld) if "647" in xx][0]
imgs = [os.path.join(cellfld, xx) for xx in os.listdir(cellfld) if "tif" in xx]
ydim,xdim = tifffile.imread(imgs[0]).shape
zdim = len(imgs)

#set dest
if not os.path.exists(os.path.join(dst, brain)): os.mkdir(os.path.join(dst, brain))
cell_dst = os.path.join(dst, brain+"/cells")
if not os.path.exists(cell_dst): os.mkdir(cell_dst)

#fill map
for zpln in range(zdim):
    if zpln%10==0: print(zpln)
    cells_dfpln = cells[cells["z"]==zpln]
    zyx = np.array([cells_dfpln["z"].values,cells_dfpln["y"].values,cells_dfpln["x"].values]).T
    cell_map = np.zeros((ydim,xdim))
    for z,y,x in zyx:
        try:
            cell_map[y,x] = 1
        except Exception as e:
            print(e)
    result = gaussian_filter(cell_map.astype(float), sigma=2)
    tifffile.imsave(os.path.join(cell_dst, "cells_%04d.tif" % zpln), 
                    result.astype(bool).astype("uint8")*255)
