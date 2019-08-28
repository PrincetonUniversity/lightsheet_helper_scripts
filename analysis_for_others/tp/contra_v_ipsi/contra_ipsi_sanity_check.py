#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 09:49:42 2019

@author: wanglab
"""

import numpy as np, tifffile as tif, cv2, matplotlib.pyplot as plt
from skimage.morphology import ball

#first, make a map of cells
pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/pma_to_aba/20180410_jg52_bl6_lob7_05/posttransformed_zyx_voxels.npy"
converted_points = np.load(pth)

zyx = np.asarray([(int(xx[0]), int(xx[1]), int(xx[2])) for xx in converted_points]) #cells are counted in horizontal volumes
   
#read vol
vol = "/jukebox/wang/zahra/h129_contra_vs_ipsi/reg_to_allen/20180410_jg52_bl6_lob7_05/cell_to_reg/result.tif"
img = tif.imread(vol)           

#init empty vol 
cell_map = np.zeros_like(img).astype(bool) 

for z,y,x in zyx:
    try:
        cell_map[z-1:z+1,y,x] = True #z dilation
    except Exception as e:
        print(e)
        
#apply x y dilation
r = 2
selem = ball(r)[int(r/2)]
cell_map = cell_map.astype("uint8")
cell_map = np.asarray([cv2.dilate(cell_map[i], selem, iterations = 1) for i in range(cell_map.shape[0])])

merged = np.stack([img, cell_map, np.zeros_like(cell_map)], -1) #rgb image you can open up in fiji; volume = red; cells = green

tif.imsave("/home/wanglab/Desktop/20180410_jg52_bl6_lob7_05_test.tif", merged)
#%%

z = [xx[0] for xx in arr]
y = [xx[1] for xx in arr]
x = [xx[2] for xx in arr]

#%%
atl_pth = "/home/wanglab/mounts/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"

#cut annotation file in middle
ann = tif.imread(atl_pth)
plt.imshow(ann[300])
#make horizontal
ann_h = np.transpose(ann, [2, 1, 0])
plt.imshow(ann_h[120])
z,y,x = ann_h.shape
ann_h_left = ann_h[:, :, :int(x/2)] #cut in the middle in x
ann_h_right = ann_h[:, :, int(x/2):]
plt.imshow(ann_h_left[120])

#%%
ann_sag_left = np.transpose(ann_h_left, [2, 1, 0])
ann_sag_right = np.transpose(ann_h_right, [2, 1, 0])

plt.imshow(ann_sag_right[0])

tif.imsave("/home/wanglab/Desktop/ann_sag_left.tif", ann_sag_left)