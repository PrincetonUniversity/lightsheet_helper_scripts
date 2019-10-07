#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 15:25:38 2019

@author: wanglab
"""

import os, SimpleITK as sitk, numpy as np, tifffile

pth = "/home/wanglab/mounts/LightSheetTransfer/microscope_tests/20190905_19_17_27_Princeton-4x_tiffs/Ex_561/ts_out/RES(1366x1799x375)/127630/127630_103390"


imgs = [os.path.join(pth, xx) for xx in os.listdir(pth)]; imgs.sort()

print(len(imgs))

arr = np.array([sitk.GetArrayFromImage(sitk.ReadImage(im)) for im in imgs])

tifffile.imsave("/home/wanglab/mounts/LightSheetTransfer/microscope_tests/20190905_19_17_27_Princeton-4x_tiffs/stitched_4x_555_3xdownsampled_stack_all_planes.tif", arr)