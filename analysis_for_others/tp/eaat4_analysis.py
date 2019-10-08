#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:45:02 2019

@author: wanglab
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np
from skimage.exposure import equalize_adapthist, rescale_intensity, adjust_gamma

img = tif.imread("/home/wanglab/mounts/wang/pisano/tracing_output/eaat4/an03_eaat4_031919/an3_eaat4_031919_1d3x_647_017na_1hfds_z10um_150msec_resized_ch00.tif")

plt.imshow(img[300]*10)

cb = img[:, 600:, :]
cb[cb < 200] = 0
plt.imshow(cb[300]*10)

z,y,x = img.shape

clh = []
for i in range(300):    
    gm = adjust_gamma(cb[i], gamma = 0.5)
    clh.append(equalize_adapthist(gm, clip_limit = 0.002, nbins = 65535))

tif.imsave("/home/wanglab/Desktop/clahe.tif", np.array(clh).astype("float32"))

plt.imshow(gm)
plt.imshow(clh)
