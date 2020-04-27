#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 18:24:12 2020

@author: wanglab
"""

import os, numpy as np, sys, matplotlib as mpl, SimpleITK as sitk
import pandas as pd, matplotlib.pyplot as plt
sys.path.append("/jukebox/wang/zahra/python/BrainPipe")
from skimage.external import tifffile
from tools.registration.transform import count_structure_lister, transformed_pnts_to_allen_helper_func
plt.ion()

src = "/jukebox/wang/Jess/lightsheet_output/202002_cfos/injection/pooled_analysis"
dst = "/jukebox/wang/Jess/lightsheet_output/202002_cfos/injection/pooled_analysis_unilateral"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
id_tb_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif",
crop = "[:, 450:, :]"
reor = ("2","0","1")
allen_id_table=pd.read_excel(id_tb_pth)
ann = sitk.GetArrayFromImage(sitk.ReadImage(ann_pth))
if crop: ann = eval("ann{}".format(crop))

#split annotation file into left and right
ann_left = np.zeros_like(ann)
zr,yr,xr = ann.shape
ann_left[:int(zr/2),:,:] = ann[:int(zr/2),:,:]
ann_right = np.zeros_like(ann)
ann_right[int(zr/2):,:,:] = ann[int(zr/2):,:,:]

plt.imshow(ann_right[300])
plt.imshow(ann_left[400])
#%%
side = "right"
#grab inj vols
nonzeros = []
vols = [os.path.join(src, xx) for xx in os.listdir(src) 
        if xx[-7:] == "inj.tif" and "crus1" in xx]
for i,vol in enumerate(vols):
    print("\n*******"+os.path.basename(vol)+"*******\n")
    nz = np.nonzero(tifffile.imread(vol))
    nonzeros.append(list(zip(*nz))) #<-for pooled image
    pos = transformed_pnts_to_allen_helper_func(np.asarray(list(zip(*[nz[2], 
                  nz[1], nz[0]]))), ann_right)
    tdf = count_structure_lister(allen_id_table, *pos)
    if i == 0: 
        df = tdf.copy()
        countcol = "count" if "count" in df.columns else "cell_count"
        df.drop([countcol], axis=1, inplace=True)
    df[os.path.basename(vol[:-7])] = tdf[countcol]

df.to_csv(os.path.join(dst,"voxel_counts_%s.csv" % side))
print("\n\nCSV file of cell counts, saved as {}\n\n\n".format(os.path.join(dst,
  "voxel_counts_%s.csv" % side)))  

side = "left"
#grab inj vols
nonzeros = []
vols = [os.path.join(src, xx) for xx in os.listdir(src) if xx[-7:] == "inj.tif"]
for i,vol in enumerate(vols):
    print("\n*******"+os.path.basename(vol)+"*******\n")
    nz = np.nonzero(tifffile.imread(vol))
    nonzeros.append(list(zip(*nz))) #<-for pooled image
    pos = transformed_pnts_to_allen_helper_func(np.asarray(list(zip(*[nz[2], 
                  nz[1], nz[0]]))), ann_left)
    tdf = count_structure_lister(allen_id_table, *pos)
    if i == 0: 
        df = tdf.copy()
        countcol = "count" if "count" in df.columns else "cell_count"
        df.drop([countcol], axis=1, inplace=True)
    df[os.path.basename(vol[:-7])] = tdf[countcol]

df.to_csv(os.path.join(dst,"voxel_counts_%s.csv" % side))
print("\n\nCSV file of cell counts, saved as {}\n\n\n".format(os.path.join(dst,
  "voxel_counts_%s.csv" % side)))      