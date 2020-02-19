#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 11:18:24 2020

@author: wanglab
"""

import tifffile, numpy as np, pandas as pd, cv2
from skimage.morphology import ball


img_pth = "/home/wanglab/wang/zahra/tracing_projects/mapping_paper/figure_data/jg51_647_z90-111_img.tif"

img = tifffile.imread(img_pth)

df = pd.read_csv("/home/wanglab/wang/zahra/tracing_projects/mapping_paper/figure_data/jg51_z90-111_cell_centers.csv")

z,y,x = df["z"].values, df["y"].values, df["x"].values

#make cnn map
cnn = np.zeros_like(img)

#map to original image
for i, zi in enumerate(z):
    cnn[zi, y[i], x[i]] = 255

#apply x y dilation
r = 4
selem = ball(r)[int(r/2)]
cnn = cnn.astype("uint8")
cnn = np.asarray([cv2.dilate(cnn[i], selem, iterations = 1) for i in range(cnn.shape[0])])

merged = np.stack([img, cnn, np.zeros_like(img)], -1)

tifffile.imsave("/home/wanglab/wang/zahra/tracing_projects/mapping_paper/figure_data/jg51_647_z90-111_cell_center_overlay.tif", merged)