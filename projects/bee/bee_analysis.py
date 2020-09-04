#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 09:47:38 2020

@author: wanglab
"""

import tifffile

#read original template
pth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_rotate.tif"
img = tifffile.imread(pth)
#crop in Z
img_mod = img[:190]
dst = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_rotate_croppedZ.tif"
tifffile.imsave(dst, img_mod)

#do same for annotation
pth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment.tif"
img = tifffile.imread(pth)
#crop in Z
img_mod = img[:190]
dst = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/template/Bombus45_2.575umstep_segment_croppedZ.tif"
tifffile.imsave(dst, img_mod)

#read in jacobian
jcpth = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/brain_to_template/Grp16_2.575_elastix/spatialJacobian.tif"
jac = tifffile.imread(jcpth)