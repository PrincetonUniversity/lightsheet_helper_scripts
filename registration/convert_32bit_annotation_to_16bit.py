#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:39:31 2019

@author: wanglab
"""

import tifffile as tif, SimpleITK as sitk, numpy as np, pandas as pd

#allen atlas
ann_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans.tif"
ann = sitk.GetArrayFromImage(sitk.ReadImage(ann_pth))

#read all unqiue structures in annotation volume
atl_ids = np.unique(ann)
atl_ids = atl_ids.astype(int)

id_table_pth = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
id_df = pd.read_excel(id_table_pth, index_col = None)
#get all ids in look up table
atl_ids_frm_tbl = id_df["id"].values.astype(int)
#find ids from the annotation volume that are also in look up table
#(should be all of them)
atl_ids_in_ann_n_tbl = [xx for xx in atl_ids if xx in id_df["id"].values]
#find names of structures that are both in annotation volume and look up table
ids_in_ann_n_tbl = id_df[id_df["id"].isin(atl_ids_in_ann_n_tbl)]
structs_in_ann_n_tbl = ids_in_ann_n_tbl["name"]

print("\n\n%d unique structures in annotation volume \n%d structures listed in look up table\
      \n%d structures in both annotation file and look up table\n" % (len(atl_ids), len(atl_ids_frm_tbl), len(structs_in_ann_n_tbl)))

#find annotation structures in annotation volume but not look up table (aka dont have a mapping)
atl_ids_in_ann_but_not_in_tbl = [xx for xx in atl_ids if not xx in id_df["id"].values] #unnamed values...

problem_ids = [xx for xx in atl_ids if xx > 65535] #range of 16 bit
#find problematic structure name, to double check things later
problem_names = [id_df.loc[id_df["id"]==xx, "name"] for xx in problem_ids]
#find other ids for these problem id regions!!!
#find highest number of ids besides these out of range ones
max_id = max([xx for xx in atl_ids if xx not in problem_ids])
#set the alterative +1 after this...
max_id += 1

ann_16bit = ann.copy()
for idx in problem_ids:
    ann_16bit[ann_16bit == idx] = max_id
    #change the identifying id in the table also
    id_df.loc[id_df["id"] == idx, "id"] = max_id
    print(id_df.loc[id_df["id"] == idx, "name"].values)
    max_id += 1

#check
new_ids = np.unique(ann_16bit)
new_ids_in_tbl = id_df["id"].values
new_ids_in_ann_n_tbl = [xx for xx in new_ids if xx in new_ids_in_tbl]

ann_16bit_converted = ann_16bit.astype("uint16")

tif.imsave("/jukebox/LightSheetTransfer/atlas/allen_atlas/annotation_2017_25um_sagittal_forDVscans_16bit.tif", ann_16bit_converted)
id_df.to_excel("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts_16bit.xlsx", index = None)

#%%
#pma
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"

ann = sitk.GetArrayFromImage(sitk.ReadImage(ann_pth))

assert ann.dtype == "float32"

atl_ids = np.unique(ann)
atl_ids = atl_ids.astype(int)

id_table_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"

id_df = pd.read_excel(id_table_pth, index_col = None)
atl_ids_frm_tbl = id_df["id"].values.astype(int)

atl_ids_in_ann_n_tbl = [xx for xx in atl_ids if xx in id_df["id"].values]

ids_in_ann_n_tbl = id_df[id_df["id"].isin(atl_ids_in_ann_n_tbl)]
structs_in_ann_n_tbl = ids_in_ann_n_tbl["name"]

print("\n\n%d unique structures in annotation volume \n%d structures listed in look up table\
      \n%d structures in both annotation file and look up table\n" % (len(atl_ids), len(atl_ids_frm_tbl), len(structs_in_ann_n_tbl)))

atl_ids_in_ann_but_not_in_tbl = [xx for xx in atl_ids if not xx in id_df["id"].values] #unnamed values...

problem_ids = [xx for xx in atl_ids if xx > 65535] #range of 16 bit

#find other ids for these problem id regions!!!
#find highest number of ids besides these out of range ones
max_id = max([xx for xx in atl_ids if xx not in problem_ids])
#set the alterative +1 after this...
max_id += 1

ann_16bit = ann.copy()
for idx in problem_ids:
    ann_16bit[ann_16bit == idx] = max_id
    #change the identifying id in the table also
    id_df.loc[id_df["id"] == idx, "id"] = max_id
    max_id += 1
    
#check
new_ids = np.unique(ann_16bit)
new_ids_in_tbl = id_df["id"].values
new_ids_in_ann_n_tbl = [xx for xx in new_ids if xx in new_ids_in_tbl]

ann_16bit_converted = ann_16bit.astype("uint16")

id_df.to_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts_16bit.xlsx", index = None)

tif.imsave("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif", ann_16bit_converted)