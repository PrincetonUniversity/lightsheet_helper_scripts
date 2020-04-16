#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:31:59 2019

@author: wanglab
"""

import os, subprocess as sp, tifffile, numpy as np, shutil, matplotlib.pyplot as plt, matplotlib as mpl
from tools.analysis.analyze_injection_inverse_transform import pool_injections_inversetransform
from tools.utils.io import makedir, load_kwargs, listdirfull
from tools.imageprocessing.orientation import fix_orientation

data = "/jukebox/wang/pisano/tracing_output/eaat4"
src = "/jukebox/wang/zahra/eaat4_screening/201910_analysis/transformed_volumes"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
dst = "/jukebox/wang/zahra/eaat4_screening/201910_analysis/merged_volumes"; makedir(dst)

imgs = listdirfull(src, "trnsfm2atl"); imgs.sort()

sites = np.array([tifffile.imread(xx) for xx in imgs]) #the y-axis cutoff for visualization

ann_raw = tifffile.imread(ann_pth) #apparent cutoff
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)

#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lingula (I)": 912,
        "Lobule II": 976,
        "Lobule III": 984,
        "Lobule IV-V": 1091,
        "Lobule VIa": 936,
        "Lobule VIb": 1134,
        "Lobule VII": 944,
        "Lobule VIII": 951,
        "Lobule IX": 957, #uvula IX
        "Lobule X": 968, #nodulus X
        "Simplex lobule": 1007, #simplex
        "Crus 1": 1056, #crus 1
        "Crus 2": 1064, #crus 2
        "Paramedian lobule": 1025, #paramedian lob
        "Copula pyramidis": 1033, #copula pyramidis
        "Paraflocculus": 1041,
        "Flocculus": 1049
        }
ak = np.array([k for k,v in iids.items()])

#find only cerbellum
cb_roi = np.zeros_like(ann_raw)
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    cb_roi[z,y,x] = 1 #make a mask of structure in annotation space

cb_roi = cb_roi.astype(bool) #add mask of structure to dictionary

#mask all the regions that are not cerebellum in the segmentation
for site in sites:
    site[~cb_roi] = 0
    
#%%
#visualization
atl = fix_orientation(tifffile.imread(atl_pth)[:, 450:, :], ("2", "0", "1"))
#masked_sites = np.array([fix_orientation(site[:, 450:, :], ("2", "0", "1")) for site in sites])

for i,site in enumerate(sites):
    masked_site = fix_orientation(site[:, 450:, :], ("2", "0", "1"))
    kwargs = load_kwargs(os.path.join(data, os.path.basename(imgs[i])[:-15]))
    cellvol = [vol for vol in kwargs["volumes"] if vol.ch_type == "cellch" or vol.ch_type == "injch"][0]
    regvol =  fix_orientation(tifffile.imread(cellvol.ch_to_reg_to_atlas)[:, 450:, :], ("2", "0", "1"))
    tifffile.imsave(os.path.join(dst, os.path.basename(imgs[i])[:-15]+"maxproj_.tif"), np.max(regvol, axis = 0))
#    merged = np.stack([regvol, masked_site, np.zeros_like(atl)], -1)
#    tifffile.imsave(os.path.join(dst, os.path.basename(imgs[i])), merged.astype("uint16"))
#
#my_cmap = eval("plt.cm.{}(np.arange(plt.cm.RdBu.N))".format("viridis"))
#my_cmap[:1,:4] = 0.0  
#my_cmap = mpl.colors.ListedColormap(my_cmap)
#my_cmap.set_under("w")
#plt.figure()
#plt.imshow(np.max(atl, axis=0), cmap="gray")
##plt.imshow(np.max(np.sum(masked_sites, axis=0), axis = 0), alpha=0.90, cmap=my_cmap); plt.colorbar(); plt.axis("off")
#plt.imshow(np.max(masked_site, axis = 0), alpha=0.90, cmap=my_cmap); plt.colorbar(); plt.axis("off")
#
#plt.savefig(os.path.join(os.path.dirname(src), "heatmap.pdf"), dpi = 300, transparent = True);
#plt.close()