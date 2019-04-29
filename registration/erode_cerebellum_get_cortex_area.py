#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 11:33:11 2019

@author: wanglab
"""

import tifffile
import numpy as np
import matplotlib.pyplot as plt

atl = tifffile.imread("/home/wanglab/mounts/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif")

zyx_scale = (20, 20, 20)
microns_to_erode = 200

#NOTE THIS ESSENTIALLY SCALES PIXEL SPACE*****
from scipy.ndimage.morphology import distance_transform_edt

distance_space_inside = distance_transform_edt(atl.astype('bool'), sampling=zyx_scale)*-1 #INSIDE
distance_space_inside = np.abs(distance_space_inside)
mask = np.copy(distance_space_inside)
mask[distance_space_inside<=microns_to_erode] = 0

#zero out edges
eann = np.copy(atl)
eann[mask==0]=0

tifffile.imsave("/home/wanglab/Desktop/annotation_sagittal_atlas_20um_iso_200um_erosion.tif", eann)
