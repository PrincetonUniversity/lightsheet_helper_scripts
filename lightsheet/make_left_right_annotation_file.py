#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 10:57:36 2019

@author: wanglab
"""

import tifffile as tif, numpy as np, matplotlib.pyplot as plt, os

dst = "/jukebox/wang/zahra/h129_contra_vs_ipsi/"

left_atl_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/horizontal_atlas_20um_iso_left_dafina_annotation.tif"

right_atl_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/horizontal_atlas_20um_iso_right_dafina_annotation.tif"

l_atl = tif.imread(left_atl_pth)

r_atl = tif.imread(right_atl_pth)

ann_pth_thal = "/jukebox/wang/pisano/Python/atlas/stepwise_erosion/annotation_sagittal_atlas_20um_iso_60um_edge_erosion_160um_ventricular_erosion.tif"

ann_pth_nc = "/jukebox/wang/pisano/Python/atlas/stepwise_erosion/annotation_sagittal_atlas_20um_iso_60um_edge_erosion_80um_ventricular_erosion.tif"

thal_ann = tif.imread(ann_pth_thal)
nc_ann = tif.imread(ann_pth_nc)

#need horizontal
thal_ann_h = np.transpose(thal_ann, [2,1,0])
nc_ann_h = np.transpose(nc_ann, [2,1,0])

plt.imshow(thal_ann_h[120])
plt.imshow(l_atl[120])

l_mask = np.where(l_atl > 0, False, True)
thal_ann_h_l = thal_ann_h.copy()
nc_ann_h_l = nc_ann_h.copy()
thal_ann_h_l[l_mask] = 0
nc_ann_h_l[l_mask] = 0
plt.imshow(thal_ann_h_l[120])

r_mask = np.where(r_atl > 0, False, True)
thal_ann_h_r = thal_ann_h.copy()
nc_ann_h_r = nc_ann_h.copy()
thal_ann_h_r[r_mask] = 0
nc_ann_h_r[r_mask] = 0
plt.imshow(thal_ann_h_r[120])

tif.imsave(os.path.join(dst, "horizontal_ann_20um_iso_left_dafina_annotation_60um_edge_erosion_160um_ventricular_erosion.tif"), thal_ann_h_l)
tif.imsave(os.path.join(dst, "horizontal_ann_20um_iso_right_dafina_annotation_60um_edge_erosion_160um_ventricular_erosion.tif"), thal_ann_h_r)
tif.imsave(os.path.join(dst, "horizontal_ann_20um_iso_right_dafina_annotation_60um_edge_erosion_80um_ventricular_erosion.tif"), nc_ann_h_r)
tif.imsave(os.path.join(dst, "horizontal_ann_20um_iso_left_dafina_annotation_60um_edge_erosion_80um_ventricular_erosion.tif"), nc_ann_h_l)