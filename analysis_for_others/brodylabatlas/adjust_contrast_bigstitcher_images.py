#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 19:43:58 2019

@author: wanglab
"""

import tifffile, matplotlib.pyplot as plt

pth = "/home/wanglab/Desktop/fused_with_content_based_blending.tif"

img = tifffile.imread(pth)
img2 = np.where(img < 150, img, img*10)
plt.imshow(img2[15])

tifffile.imsave("/home/wanglab/Desktop/fused_with_content_based_blending_int_adj.tif", img2)
