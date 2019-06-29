#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 18:18:32 2019

@author: wanglab
"""

import tifffile, numpy as np
from scipy.ndimage import zoom

img_pth = "/home/wanglab/mounts/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"

vol = tifffile.imread(img_pth)
#org dims
z,y,x = vol.shape

#new dims 
zr, yr, xr = (75,137,75)

#resize 3d 
resz = zoom(vol.astype("uint16"), (zr/z, yr/y, xr/x), order = 3)

#save out 
tifffile.imsave("/jukebox/wang/zahra/python/lightsheet_py3/supp_files/ann_downsampled.tif", resz, compress = 6)