#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:27:56 2019

@author: wanglab
"""

import tifffile, matplotlib.pyplot as plt, numpy as np
from scipy.ndimage.filters import gaussian_filter as gfilt
from scipy.ndimage import label
import matplotlib.colors
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "Red"])
from tools.imageprocessing.orientation import fix_orientation
from tools.conv_net.input.read_roi import read_roi_zip
from tools.registration.transform import transformed_pnts_to_allen_helper_func
from tools.registration.register import count_structure_lister
from tools.utils.io import load_kwargs
impth = "/home/wanglab/mounts/wang/oostland/lightsheet/m26/marlies_m26_1d3x_488_555_008na_1hfds_z5um_150msec_resized_ch01_resampledforelastix.tif"

im = tifffile.imread(impth)
imcor = fix_orientation(im, axes = ("2", "0", "1")) #coronal

plt.imshow(np.max(im[:, 600:, :]*10, axis=0))

imcor = imcor[600:, :, :] #crop

#try analyze injection route
thresh=4; filter_kernel=(3,3,3)
filtered = gfilt(imcor, filter_kernel)
thresholded = filtered > filtered.mean() + thresh*filtered.std() 
labelled,nlab = label(thresholded)
num_sites_to_keep = 20

if nlab == 0:
    raise Exception("Site not detected, try a lower threshold?")
elif nlab == 1:
    arr = labelled.astype(bool)
elif num_sites_to_keep == 1:
    sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
    arr = labelled == np.argmax(sizes)+1
else:
    sizes = [np.sum(labelled==i) for i in range(1,nlab+1)]
    vals = [i+1 for i in np.argsort(sizes)[-num_sites_to_keep:][::-1]]
    arr = np.in1d(labelled, vals).reshape(labelled.shape)

#plot max proj
plt.imshow(np.max(imcor*20, axis=0), "gist_yarg")
plt.imshow(np.max(arr, axis=0), cmap, alpha = 0.5)
tifffile.imsave('/home/wanglab/Desktop/result.tif', arr)
