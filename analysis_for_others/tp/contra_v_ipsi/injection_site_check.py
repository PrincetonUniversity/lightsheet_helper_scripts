#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 17:29:12 2019

@author: wanglab
"""

import tifffile as tif, matplotlib.pyplot as plt, numpy as np, os, matplotlib as mpl

src = "/home/wanglab/mounts/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites"
atl = "/home/wanglab/mounts/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
dst = "/home/wanglab/Desktop"

#get coronal volumes
imgs = np.array([np.rot90(np.transpose(tif.imread(os.path.join(src, xx)), [1, 0, 2]), axes = (2,1)) for xx in os.listdir(src)])
atl = tif.imread(atl)
atl = atl[:, 423:, :]

#make coronal
atl_cor = np.rot90(np.transpose(atl, [1, 0, 2]), axes = (2,1))

cmap = "plasma"
my_cmap = eval("plt.cm.{}(np.arange(plt.cm.RdBu.N))".format(cmap))
my_cmap[:1,:4] = 0.0  
my_cmap = mpl.colors.ListedColormap(my_cmap)
my_cmap.set_under("w")
plt.figure()
plt.imshow(np.max(atl_cor, axis=0), cmap="gray")
plt.imshow(np.max(np.max(imgs, axis=0), axis=0), alpha=0.99, cmap=my_cmap); plt.colorbar(); plt.axis("off")
plt.savefig(os.path.join(dst,"heatmap.pdf"), dpi=300, transparent=True);
plt.close()