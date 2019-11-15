#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:54:08 2019

@author: wanglab
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, pickle as pckl, matplotlib as mpl, statsmodels.api as sm, itertools
import ast
from skimage.external import tifffile
from scipy.ndimage.measurements import center_of_mass

from tools.registration.register import transformed_pnts_to_allen_helper_func, count_structure_lister, change_transform_parameter_initial_transform
from tools.registration.transform_list_of_points import create_text_file_for_elastix, modify_transform_files, point_transformix, unpack_pnts
from tools.utils.io import makedir
from tools.analysis.network_analysis import make_structure_objects

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#set paths, variables, etc.
dst = "/jukebox/wang/zahra/prv/"
pma_ann_pth = os.path.join(dst, "pma_annotation_sagittal_atlas_20um_iso_60um_edge_160um_vntric_erosion.tif")
pma_ann = tifffile.imread(pma_ann_pth)

#collect 
#brains should be in this order as they were saved in this order for inj analysis
brains = ["20180205_jg_bl6f_prv_02", "20180205_jg_bl6f_prv_03", "20180205_jg_bl6f_prv_04", "20180215_jg_bl6f_prv_05", "20180215_jg_bl6f_prv_06",
       "20180215_jg_bl6f_prv_09", "20180305_jg_bl6f_prv_12", "20180305_jg_bl6f_prv_15", "20180312_jg_bl6f_prv_17", "20180326_jg_bl6f_prv_37",
       "20180313_jg_bl6f_prv_21", "20180313_jg_bl6f_prv_23", "20180313_jg_bl6f_prv_24", "20180305_jg_bl6f_prv_11", "20180313_jg_bl6f_prv_25",
       "20180322_jg_bl6f_prv_27", "20180322_jg_bl6f_prv_28", "20180323_jg_bl6f_prv_30", "20180326_jg_bl6f_prv_33", 
       "20180326_jg_bl6f_prv_34", "20180326_jg_bl6f_prv_35"]
    
inj_src = os.path.join(dst, "prv_injection_sites")
imgs = [os.path.join(inj_src, xx+".tif") for xx in brains]

#pool brain names and L/R designation into dict
lr_dist = {}
inj_vox = {}

#get inj vol roundabout way
for img in imgs:
    brain = os.path.basename(img)
    print(brain)
    inj_vol = tifffile.imread(img)
    z,y,x = inj_vol.shape
    
    z_c,y_c,x_c = center_of_mass(inj_vol)
    #take distance from center to arbitrary "midline" (aka half of z axis)
    dist = z_c-(z/2)
    #save to dict 
    lr_dist[brain[:-4]] = dist
    inj_vox[brain[:-4]] = inj_vol
    
    if dist < 0:
        print("brain {} has a left-sided injection\n".format(brain))
    elif dist > 0:
        print("brain {} has a right-sided injection\n".format(brain))
    else:
        print("brain has an injection close to midline so not considering it rn\n")


#get injection fractions
inj_raw = np.array([inj.astype(bool) for inj in inj_vox.values()])
    
atl_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif")
ann_raw = tifffile.imread("/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso_16bit.tif")
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
        "Copula pyramidis": 1033 #copula pyramidis
        }
ak = np.asarray([k for k,v in iids.items()])

atlas_rois = {}
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    ann_blank = np.zeros_like(ann_raw)
    ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
    atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary

expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

#pooled injections
ak_pool = np.array(["Lob. I-III, IV-V", "Lob. VIa, VIb, VII", "Lob. VIII, IX, X",
                 "Simplex", "Crura", "PM, CP"])
frac_of_inj_pool = np.array([[np.sum(xx[:4]),np.sum(xx[4:7]),np.sum(xx[7:10]),xx[10],xx[11:13].sum(),np.sum(xx[13:16])] for xx in expr_all_as_frac_of_inj])
primary_pool = np.array([np.argmax(e) for e in frac_of_inj_pool])
#get n"s after pooling
primary_lob_n = np.array([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])


#make structures
#FIXME: for some reason the allen table does not work on this, is it ok to use PMA        
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=pma_ann_pth)

#set variables
lr_brains = list(lr_dist.keys())

atl_dst = os.path.join(dst, "pma_to_aba"); makedir(atl_dst)
trnsfrm_dst = os.path.join(dst, "prv_transformed_cells")
id_table = pd.read_excel(df_pth)

#%%

def transformed_cells_to_ann(fld):
    """ consolidating to one function """
    
    trnsfrmd_pnts = []
    
    for fl in fld:
        pnts_pths = os.path.join(fl, "posttransformed_zyx_voxels.npy")
        if not os.path.exists(pnts_pths): 
            pnts_pths = os.path.join(fl, "transformed_points/posttransformed_zyx_voxels.npy")
            trnsfrmd_pnts.append(np.load(pnts_pths).astype(int))
    
    return np.array(trnsfrmd_pnts)

#pma2aba_transformed = [os.path.join(atl_dst, xx) for xx in lr_brains]
trnsfrmd = [os.path.join(trnsfrm_dst, xx) for xx in lr_brains]

whl_brains = transformed_cells_to_ann(trnsfrmd)

#%%
all_coords = []
#get all unique coordinates, see how many unique there are...
for brain in whl_brains:
    all_coords.append(list(np.unique([str(tuple(xx)) for xx in brain])))
    
all_coords = list(itertools.chain.from_iterable(all_coords))
#all unique coords
all_unique_coords = np.unique(all_coords)

print(len(all_coords))
print(len(all_unique_coords))

#jitter of 4 px (80 um)
coords = {}
i = 0
for unq in all_unique_coords:
    if not ast.literal_eval(unq)[0] < 0 or ast.literal_eval(unq)[1] < 0 or ast.literal_eval(unq)[2] < 0: #not include negative coordinates
        coords[unq] = i
        i+=1