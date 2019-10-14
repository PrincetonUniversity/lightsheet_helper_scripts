#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:45:02 2019

@author: wanglab
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np, os
from skimage.exposure import equalize_adapthist
from scipy.ndimage.filters import gaussian_filter as gfilt
from scipy.ndimage import label
from tools.utils.io import makedir, load_kwargs, listdirfull, listall
from tools.registration.transform_list_of_points import modify_transform_files, point_transformix, unpack_pnts, create_text_file_for_elastix
from tools.imageprocessing.orientation import fix_orientation
from tools.registration.transform_cell_counts import change_transform_parameter_initial_transform, points_resample

src = "/home/wanglab/mounts/wang/pisano/tracing_output/eaat4/an19_eaat4_031919/full_sizedatafld/an19_eaat4_031919_1d3x_647_017na_1hfds_z10um_150msec_ch00"

imgs = [os.path.join(src, xx) for xx in os.listdir(src)]; imgs.sort()

stk = np.array([tif.imread(img) for img in imgs])[:, 1700:, :]

clh = np.array([equalize_adapthist(img, kernel_size = (50,100), clip_limit=0.05, nbins = 65535) for img in stk])

power = 3 #how many times to multiple image by itself

plt.figure()
plt.imshow(clh[300]**power)

thresh = 2
filter_kernel=(5,5,5)

filtered = gfilt(clh**power, filter_kernel)
tif.imsave("/home/wanglab/Desktop/filtered_z200-300.tif", filtered[200:300].astype("float32"))

plt.figure()
plt.imshow(filtered[300], "gist_yarg")

thresholded = filtered > filtered.mean() + thresh*filtered.std() 
labelled,nlab = label(thresholded)

plt.figure()
sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
vals = [i+1 for i in np.argsort(sizes)[::-1]]
arr = np.in1d(labelled, vals).reshape(labelled.shape)

tif.imsave("/home/wanglab/Desktop/seg.tif", np.in1d(labelled, vals).reshape(labelled.shape).astype("uint16"))

#%%
#collect nonzero pixels

import cv2
from scipy.ndimage import zoom

arr = tif.imread("/home/wanglab/Desktop/seg.tif")
im = tif.imread(imgs[0])
arr_fullsz = np.zeros((len(imgs), im.shape[0], im.shape[1]))
print(arr_fullsz.shape)
print(arr.shape)

arr_fullsz[:, 1700:, :] = arr

plt.imshow(arr_fullsz[300])

zf,yf,xf = 1,(1/3),(1/3)
dim = (int(im.shape[1]*xf), int(im.shape[0]*yf))
arr_resz = np.array([cv2.resize(img, dim) for img in arr_fullsz])
print(arr_resz.shape)

#to sagittal
arr_resz_sag = fix_orientation(arr_resz, axes = ("2", "1", "0"))

tif.imsave("/home/wanglab/Desktop/seg_resized.tif", arr_resz_sag.astype("uint16")) 

#%%
arr_resz_sag = tif.imread("/home/wanglab/Desktop/seg_resized.tif")
#now downsample again
resmpl = tif.imread("/home/wanglab/mounts/wang/pisano/tracing_output/eaat4/an19_eaat4_031919/an19_eaat4_031919_1d3x_647_017na_1hfds_z10um_150msec_resized_ch00_resampledforelastix.tif")

zrsm, yrsm, xrsm = resmpl.shape
zresz, yresz, xresz = arr_resz_sag.shape

arr_resmpl_sag = zoom(arr_resz_sag, (zrsm/zresz, yrsm/yresz, xrsm/xresz), order = 1)

tif.imsave("/home/wanglab/Desktop/seg_downsampled.tif", arr_resmpl_sag.astype("uint16")) 
