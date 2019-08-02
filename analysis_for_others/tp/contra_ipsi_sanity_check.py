#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 09:49:42 2019

@author: wanglab
"""

import numpy as np, tifffile as tif, cv2, matplotlib.pyplot as plt
from skimage.morphology import ball

#first, make a map of cells
zyx = np.asarray([(int(xx[0]), int(xx[1]), int(xx[2])) for xx in arr]) #cells are counted in horizontal volumes
   
vol = "/jukebox/wang/pisano/tracing_output/antero_4x/20180608_jg72/elastix/20180608_bl6_h129_jg72_4x_647_008na_1hfds_z7d5um_300msec_10povlp_resized_ch00/result.tif"

img = tif.imread(vol)           
 
cell_map = np.zeros_like(img).astype(bool) 

for z,y,x in zyx:
    try:
        cell_map[z-1:z+1,y,x] = True
    except Exception as e:
        print(e)
        
#apply x y dilation
r = 2
selem = ball(r)[int(r/2)]
cell_map = cell_map.astype("uint8")
cell_map = np.asarray([cv2.dilate(cell_map[i], selem, iterations = 1) for i in range(cell_map.shape[0])])

merged = np.stack([img, cell_map, np.zeros_like(cell_map)], -1)

tif.imsave("/home/wanglab/Desktop/test.tif", merged)
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