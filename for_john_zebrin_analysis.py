#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 16:25:05 2019

@author: wanglab
"""

from skimage.external import tifffile
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np

img_pth = "/home/wanglab/mounts/wang/pisano/tracing_output/eaat4/an11_eaat4_031919/elastix/an11_eaat4_031919_1d3x_647_017na_1hfds_z10um_150msec_resized_ch00/result.tif"

#read the image
img = tifffile.imread(img_pth)

#look at one slice
plt.imshow(img[200])

#crop the image
cb = img[:, 423:, :]

#show cropped image
plt.imshow(cb[400])

#reslice to posterior-anterior coronal view
cb_coronal = np.rot90(np.transpose(cb, [1, 0, 2]), axes = (2,1))

#show rotated image
plt.imshow(cb_coronal[0])

#get max projection
max_cb = np.max(cb_coronal, axis = 0)

#show max projection
plt.imshow(max_cb, cmap = "gist_yarg")

#export as tif
tifffile.imsave("/home/wanglab/Desktop/test.tif", max_cb)
