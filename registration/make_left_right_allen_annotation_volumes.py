#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:08:16 2019

@author: wanglab
"""

import SimpleITK as sitk,  numpy as np, matplotlib.pyplot as plt, os, tifffile

pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans.nrrd"

ann = sitk.GetArrayFromImage(sitk.ReadImage(pth))

plt.imshow(ann[300])

z,y,x = ann.shape
#make sure each halves are same dimension as original ann
ann_left = np.zeros_like(ann)
ann_left[:int(z/2), :, :] = ann[:int(z/2), :, :] #cut in the middle in x
ann_right = np.zeros_like(ann)
ann_right[int(z/2):, :, :] = ann[int(z/2):, :, :]
plt.imshow(ann_left[120])

tifffile.imsave(os.path.join(os.path.dirname(pth), "annotation_2017_25um_sagittal_forDVscans_left_side_only.tif"), ann_left)
