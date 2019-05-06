#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 14:27:51 2019

@author: wanglab
"""

import tifffile, numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches as mpatches

#custom
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_200um_edges_only.tif"

atl_raw = tifffile.imread(atl_pth)[:, 423:, :]
ann_raw = tifffile.imread(ann_pth)[:, 423:, :] #apparent cutoff
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)

#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lobule IV-V": 1091,
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
        "Copula pyramidis": 1033 #copula pyramidis
        }
ak = np.asarray([k for k,v in iids.items()])

atlas_rois = {}
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    ann_blank = np.zeros_like(ann_raw)
    ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
    atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary

#%%
#vermis
fig, ax = plt.subplots()    
plt.axis('off')
ax = plt.imshow(np.max(atl_raw, axis = 0), cmap = "binary")
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "red"])
ax = plt.imshow(np.max(atlas_rois["Lobule IV-V"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "lime"])
ax = plt.imshow(np.max(atlas_rois["Lobule VIa"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "navy"])
ax = plt.imshow(np.max(atlas_rois["Lobule VIb"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "gold"])
ax = plt.imshow(np.max(atlas_rois["Lobule VII"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "coral"])
ax = plt.imshow(np.max(atlas_rois["Lobule VIII"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "indigo"])
ax = plt.imshow(np.max(atlas_rois["Lobule IX"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "fuchsia"])
ax = plt.imshow(np.max(atlas_rois["Lobule X"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "chocolate"])
ax = plt.imshow(np.max(atlas_rois["Copula pyramidis"].astype(int), axis = 0), cmap = cmap, alpha = 0.4)

red_patch = mpatches.Patch(color="red", label="Lobule IV-V")
green_patch = mpatches.Patch(color="lime", label="Lobule VIa")
navy_patch = mpatches.Patch(color="navy", label="Lobule VIb")
gold_patch = mpatches.Patch(color="gold", label="Lobule VII")
coral_patch = mpatches.Patch(color="coral", label="Lobule VIII")
azure_patch = mpatches.Patch(color="indigo", label="Lobule IX")
fuchsia_patch = mpatches.Patch(color="fuchsia", label="Lobule X")
maroon_patch = mpatches.Patch(color="chocolate", label="Copula pyramidis")

plt.legend(handles=[red_patch, green_patch, navy_patch, gold_patch, 
                    coral_patch, azure_patch, fuchsia_patch, maroon_patch], bbox_to_anchor=(.8, 1), loc=2, borderaxespad=0.)

plt.savefig("/home/wanglab/Desktop/test.svg")

#%%
#hemisphere
fig, ax = plt.subplots()    
plt.axis('off')
ax = plt.imshow(atl_raw[150], cmap = "binary")
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "aquamarine"])
ax = plt.imshow(atlas_rois["Simplex lobule"].astype(int)[150], cmap = cmap, alpha = 0.5)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "chartreuse"])
ax = plt.imshow(atlas_rois["Crus 1"].astype(int)[150], cmap = cmap, alpha = 0.5)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "salmon"])
ax = plt.imshow(atlas_rois["Crus 2"].astype(int)[150], cmap = cmap, alpha = 0.5)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "violet"])
ax = plt.imshow(atlas_rois["Paramedian lobule"].astype(int)[150], cmap = cmap, alpha = 0.5)

aquamarine_patch = mpatches.Patch(color="aquamarine", label="Simplex lobule")
chartreuse_patch = mpatches.Patch(color="chartreuse", label="Crus 1")
salmon_patch = mpatches.Patch(color="salmon", label="Crus 2")
violet_patch = mpatches.Patch(color="violet", label="Paramedian lobule")

plt.legend(handles=[aquamarine_patch, chartreuse_patch, 
                    salmon_patch, violet_patch], bbox_to_anchor=(.8, 1), loc=2, borderaxespad=0.)

plt.savefig("/home/wanglab/Desktop/test.svg")
