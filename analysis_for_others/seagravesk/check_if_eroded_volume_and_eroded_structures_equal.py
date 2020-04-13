#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 17:06:17 2020

@author: wanglab
"""

import os, pandas as pd, tifffile

csvpth = "/jukebox/wang/zahra/kelly_cell_detection_analysis/structures_zeroed_after_80um_erosion_allen_annotation_2017_ids_fixed_v2.csv"
df = pd.read_csv(csvpth)

erode_pth = "/jukebox/wang/zahra/kelly_cell_detection_analysis/annotation_allen_2017_25um_sagittal_erode_80um.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans_16bit.tif"

erode = tifffile.imread(erode_pth)
ann = tifffile.imread(ann_pth)

erode_iid = np.unique(erode).astype(int)
ann_iid = np.unique(ann).astype(int)

structs_eroded = [xx for xx in ann_iid if xx not in erode_iid]
print(len(structs_eroded))

df_iids_eroded = df["id"].values
print(len(df))

missing_frm_df = [xx for xx in structs_eroded if xx not in df_iids_eroded]
print(missing_frm_df)