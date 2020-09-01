#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 09:47:38 2020

@author: wanglab
"""

import os, numpy as np, sys, time
import tifffile, SimpleITK as sitk
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from tools.utils.io import makedir, load_memmap_arr, listall, load_kwargs
from tools.registration.register import change_interpolation_order, transformix_command_line_call
from tools.registration.transform_list_of_points import modify_transform_files
from scipy.ndimage.interpolation import zoom

tr0 = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575_elastix/TransformParameters.0.txt"
tr1 = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575_elastix_2/TransformParameters.0.txt"
dst = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575_transform"
if not os.path.exists(dst): os.mkdir(dst)
transformfiles = modify_transform_files(transformfiles=[tr0,tr1], dst=dst)
[change_interpolation_order(xx,0) for xx in transformfiles]

#%%

#read jacobian determinant
imgpth1 = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/IsoYellow_2.575_jac/spatialJacobian.tif"
jac1 = tifffile.imread(imgpth1)

imgpth2 = "/home/wanglab/LightSheetData/kocher-bee/volume_analysis/Grp16_2.575_jac/spatialJacobian.tif"
jac2 = tifffile.imread(imgpth2)

