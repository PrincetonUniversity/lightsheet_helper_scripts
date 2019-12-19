#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 14:25:05 2019

@author: wanglab
"""

import sys
sys.path.append("/jukebox/wang/zahra/python/lightsheet_helper_scripts")
from lightsheet.network_analysis import make_structure_objects
from scipy.stats import spearmanr
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, os, pickle as pckl, pandas as pd
from skimage.external import tifffile
import seaborn as sns

#custom
inj_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/linear_modeling/thalamus/injection_sites"
atl_pth = "/jukebox/LightSheetTransfer/atlas/sagittal_atlas_20um_iso.tif"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
cells_regions_pth = "/jukebox/wang/zahra/h129_contra_vs_ipsi/data/thal_contra_counts_23_brains.csv"
dst = "/home/wanglab/Desktop"
#making dictionary of injection sites
injections = {}

#MAKE SURE THEY ARE IN THIS ORDER
brains = ['20170410_tp_bl6_lob6a_ml_repro_01',
#         '20160823_tp_bl6_cri_500r_02',
         '20180417_jg59_bl6_cri_03',
         '20170207_db_bl6_crii_1300r_02',
         '20160622_db_bl6_unk_01',
         '20161205_tp_bl6_sim_750r_03',
         '20180410_jg51_bl6_lob6b_04',
         '20170419_db_bl6_cri_rpv_53hr',
         '20170116_tp_bl6_lob6b_lpv_07',
         '20170411_db_bl6_crii_mid_53hr',
         '20160822_tp_bl6_crii_1500r_06',
         '20160920_tp_bl6_lob7_500r_03',
         '20170207_db_bl6_crii_rpv_01',
         '20161205_tp_bl6_sim_250r_02',
         '20161207_db_bl6_lob6a_500r_53hr',
         '20170130_tp_bl6_sim_rlat_05',
         '20170115_tp_bl6_lob6b_500r_05',
         '20170419_db_bl6_cri_mid_53hr',
         '20161207_db_bl6_lob6a_850r_53hr',
         '20160622_db_bl6_crii_52hr_01',
         '20161207_db_bl6_lob6a_50rml_53d5hr',
         '20161205_tp_bl6_lob45_1000r_01',
         '20160801_db_l7_cri_01_mid_64hr']

for pth in brains:
    print(pth)
    pth = os.path.join(inj_pth, pth+".tif.tif")
    injection = tifffile.imread(pth)
    print(injection.shape)
    injections[os.path.basename(pth)[:-8]] = injection #files have 2 .tif in the end
    
inj_raw = np.array([inj.astype(bool) for nm, inj in injections.items()])
    
atl_raw = tifffile.imread(atl_pth)
ann_raw = tifffile.imread(ann_pth)[:, 423:, :] #apparent cutoff
anns = np.unique(ann_raw).astype(int)
print(ann_raw.shape)

#annotation IDs of the cerebellum ONLY that are actually represented in annotation file
iids = {"Lobule III": 984,
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

#just cerebellar region names
ak = np.asarray([k for k,v in iids.items()])

atlas_rois = {}
for nm, iid in iids.items():
    z,y,x = np.where(ann_raw == iid) #find where structure is represented
    ann_blank = np.zeros_like(ann_raw)
    ann_blank[z,y,x] = 1 #make a mask of structure in annotation space
    atlas_rois[nm] = ann_blank.astype(bool) #add mask of structure to dictionary
    
#get fractions
expr_all_as_frac_of_lob = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / lob.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
expr_all_as_frac_of_inj = np.array([[(mouse.ravel()[lob.ravel()].astype(int)).sum() / mouse.sum() for nm, lob in atlas_rois.items()] for mouse in inj_raw])    
primary = np.array([np.argmax(e) for e in expr_all_as_frac_of_inj])
primary_as_frac_of_lob = np.array([np.argmax(e) for e in expr_all_as_frac_of_lob])
secondary = np.array([np.argsort(e)[-2] for e in expr_all_as_frac_of_inj])

print(expr_all_as_frac_of_lob[15])

#%%
#making dictionary of cells by region
cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
#pooling regions
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)


#%%
#atlas res
scale_factor = 0.020 #mm/voxel

#GET ONLY VPM + bonus thal nuclei?, VTA, AND SNC COUNTS TO COMPARE

cells_regions = pd.read_csv(cells_regions_pth)
#rename structure column
cells_regions["Structure"] = cells_regions["Unnamed: 0"]
cells_regions = cells_regions.drop(columns = ["Unnamed: 0"])
ann_df = "/jukebox/LightSheetTransfer/atlas/allen_atlas/allen_id_table_w_voxel_counts.xlsx"
scale_factor = 0.025
ann_df = pd.read_excel(ann_df).drop(columns = ["Unnamed: 0"])

#get counts for all of neocortex
sois = ["Ventral tegmental area", #vta
        "Substantia nigra, reticular part", 
        "Substantia nigra, compact part",#snc
        "Reticular nucleus of the thalamus", #thal
        "Mediodorsal nucleus of thalamus",
        "Ventral posteromedial nucleus of the thalamus"
        ]


#make new dict - for all brains
cells_pooled_regions = {} #for raw counts
vol_pooled_regions = {}

for brain in brains:    
    #make new dict - this is for EACH BRAIN
    c_pooled_regions = {}
    v_pooled_regions = {}
    
    for soi in sois:
        counts = []; vol = []
        soi = [s for s in structures if s.name==soi][0]
        
        counts.append(cells_regions.loc[cells_regions.Structure == soi.name, brain].values[0]) #store counts in this list
        vol.append(ann_df.loc[ann_df.name == soi.name, "voxels_in_structure"].values[0]) #store vols in this list
        #add to volume list from LUT
        progeny = [str(xx.name) for xx in soi.progeny]
        #now sum up progeny
        if len(progeny) > 0:
            for progen in progeny:
                counts.append(cells_regions.loc[cells_regions.Structure == progen, brain].values[0])
                vol.append(ann_df.loc[ann_df.name == progen, "voxels_in_structure"].values[0])
    
        c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
        v_pooled_regions[soi.name] = np.sum(np.asarray(vol))
    
    #add to big dict
    cells_pooled_regions[brain] = c_pooled_regions
    vol_pooled_regions[brain] = v_pooled_regions
    
#making the proper array per brain where regions are removed
cell_counts_per_brain = []

#initialise dummy var
i = []
for k,v in cells_pooled_regions.items():
    dct = cells_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    cell_counts_per_brain.append(i)
    #re-initialise for next
    i = []

#making the proper array per brain where regions are removed
vol_per_brain = []

#initialise dummy var
i = []
for k,v in vol_pooled_regions.items():
    dct = vol_pooled_regions[k]
    for j,l in dct.items():
        i.append(l)  
    vol_per_brain.append(i)
    #re-initialise for next
    i = []
    
cell_counts_per_brain = np.asarray(cell_counts_per_brain)
vol_per_brain = np.array(vol_per_brain)/2
density_per_brain = np.asarray([xx/(vol_per_brain[i]*(scale_factor**3)) for i, xx in enumerate(cell_counts_per_brain)])

#rename
short_nuclei = ["VTA", #vta
        "SNc", #snc
        "SNr", #snc
        "RTN", #thal
        "MD",
        "VPM"
        ]

#%%
#pooled injections
ak_pool = np.asarray(["Lob. III, IV-V", "Lob. VIa, VIb, VII-X", 
                 "Simplex", "Crus I", "Crus II", "PM, CP"])

#pooling injection regions
expr_all_as_frac_of_lob_pool = np.asarray([[brain[0]+brain[1], brain[2]+brain[3]+brain[4]+brain[5]+brain[6]+brain[7], 
                                            brain[8], brain[9], brain[10], brain[11]+brain[12]] for brain in expr_all_as_frac_of_lob])

expr_all_as_frac_of_inj_pool = np.asarray([[brain[0]+brain[1], brain[2]+brain[3]+brain[4]+brain[5]+brain[6]+brain[7], 
                                            brain[8], brain[9], brain[10], brain[11]+brain[12]] for brain in expr_all_as_frac_of_inj])
    
primary_pool = np.asarray([np.argmax(e) for e in expr_all_as_frac_of_inj_pool])
primary_lob_n = np.asarray([np.where(primary_pool == i)[0].shape[0] for i in np.unique(primary_pool)])

print(expr_all_as_frac_of_lob_pool.shape)

#normalise inj
total_lob_sum = np.asarray([np.sum(expr_all_as_frac_of_lob_pool[i]) for i in range(len(expr_all_as_frac_of_lob))])    
expr_all_as_frac_of_lob_pool_norm = np.asarray([expr_all_as_frac_of_lob_pool[i]/total_lob_sum[i] for i in range(len(expr_all_as_frac_of_lob))])

#%%
## CELL COUNTS
#ignoring cb topology as this can be confusing and non-specific
fig = plt.figure(figsize=(15,2))
ax = fig.add_axes([.4,.1,.5,.8])

show = cell_counts_per_brain.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 100
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", 
                  ticks=bounds, boundaries=bounds, format="%d", shrink=0.5, aspect=10)
cb.set_label("Cell counts", fontsize="x-small", labelpad=2)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < 40:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}\ninj={}".format(br, ak) for br, ak in zip(brains, ak_pool[primary_pool])], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)

ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="xx-small")

plt.savefig(os.path.join(dst,"thalvtacomp_cell_count.pdf"), bbox_inches = "tight")

#%%    
#sort density
sorted_counts = [cell_counts_per_brain[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sorted_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]

#reformat - wtf
import itertools
sorted_counts = np.asarray(list(itertools.chain.from_iterable(sorted_counts)))
sorted_brains = list(itertools.chain.from_iterable(sorted_brains))
sorted_inj = np.asarray(['Lob. III, IV-V', 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X',
 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X',
 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X', 'Simplex', 'Simplex', 'Crus I', 'Crus I', 'Crus I', 'Crus I', 
 'Crus II', 'Crus II', 'Crus II', 'Crus II', 'PM, CP', 'PM, CP', 'PM, CP', 'PM, CP'])
    
fig = plt.figure(figsize=(15,2))
ax = fig.add_axes([.4,.1,.5,.8])

vmin = 0
vmax = 100
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", 
                  ticks=bounds, boundaries=bounds, format="%d", shrink=0.5, aspect=10)
cb.set_label("Cell counts", fontsize="x-small", labelpad=2)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < 40:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}\ninj={}".format(br, ak) for br, ak in zip(sorted_brains, sorted_inj)], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)

ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="xx-small")

plt.savefig(os.path.join(dst,"thalvtacomp_cell_counts_sorted.pdf"), bbox_inches = "tight")

#%%
## DENSITY
#ignoring cb topology as this can be confusing and non-specific
fig = plt.figure(figsize=(15,2))
ax = fig.add_axes([.4,.1,.5,.8])

show = density_per_brain.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 75
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", 
                  ticks=bounds, boundaries=bounds, format="%0.1f", shrink=0.5, aspect=10)
cb.set_label("Cells/mm3", fontsize="x-small", labelpad=2)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < 30:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}\ninj={}".format(br, ak) for br, ak in zip(brains, ak_pool[primary_pool])], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)

ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="xx-small")

plt.savefig(os.path.join(dst,"thalvtacomp_density.pdf"), bbox_inches = "tight")

#%%
#sort density
sorted_counts = [density_per_brain[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]
sorted_brains = [np.asarray(brains)[np.where(primary_pool == idx)[0]] for idx in np.unique(primary_pool)]

#reformat - wtf
import itertools
sorted_counts = np.asarray(list(itertools.chain.from_iterable(sorted_counts)))
sorted_brains = list(itertools.chain.from_iterable(sorted_brains))
sorted_inj = np.asarray(['Lob. III, IV-V', 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X',
 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X',
 'Lob. VIa, VIb, VII-X', 'Lob. VIa, VIb, VII-X', 'Simplex', 'Simplex', 'Crus I', 'Crus I', 'Crus I', 'Crus I', 
 'Crus II', 'Crus II', 'Crus II', 'Crus II', 'PM, CP', 'PM, CP', 'PM, CP', 'PM, CP'])
    
fig = plt.figure(figsize=(15,2))
ax = fig.add_axes([.4,.1,.5,.8])

show = sorted_counts.T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 75
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", 
                  ticks=bounds, boundaries=bounds, format="%0.1f", shrink=0.5, aspect=10)
cb.set_label("Cells/mm3", fontsize="x-small", labelpad=2)
cb.ax.tick_params(labelsize="x-small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < 30:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="white", ha="center", va="center", fontsize="xx-small")
        else:
            ax.text(ci+.5, ri+.5, "{:0.1f}".format(col), color="k", ha="center", va="center", fontsize="xx-small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(sorted_brains))+.5)
lbls = np.asarray(brains)
ax.set_xticklabels(["{}\ninj={}".format(br, ak) for br, ak in zip(sorted_brains, sorted_inj)], rotation=45, fontsize=5, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)

ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="xx-small")

plt.savefig(os.path.join(dst,"thalvtacomp_density_sorted.pdf"), bbox_inches = "tight")
    
#%%
## mean density
fig = plt.figure(figsize=(5,2))
ax = fig.add_axes([.4,.1,.5,.8])

mean_density = np.asarray([np.mean(density_per_brain[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

show = mean_density.astype(int).T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 75
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", 
                  ticks=bounds, boundaries=bounds, format="%d", shrink=0.5, aspect=10)
cb.set_label("Cells/$mm^3$", fontsize="small", labelpad=2)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < 20:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(ak_pool, primary_lob_n)], rotation=45, fontsize=7, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)

ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="small")

plt.savefig(os.path.join(dst,"thalvtacomp_mean_density.pdf"), bbox_inches = "tight")

#%%

## mean counts
fig = plt.figure(figsize=(5,2))
ax = fig.add_axes([.4,.1,.5,.8])

mean_counts = np.asarray([np.mean(cell_counts_per_brain[np.where(primary_pool == idx)[0]], axis=0) for idx in np.unique(primary_pool)])

show = mean_counts.astype(int).T #np.flip(mean_counts, axis = 1) # NOTE abs

vmin = 0
vmax = 50
cmap = plt.cm.viridis
cmap.set_over('gold')
#colormap
# discrete colorbar details
bounds = np.linspace(vmin,vmax,6)
#bounds = np.linspace(-2,5,8)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pc = ax.pcolor(show, cmap=cmap, vmin=vmin, vmax=vmax)#, norm=norm)
cb = plt.colorbar(pc, ax=ax, cmap=cmap, norm=norm, spacing="proportional", 
                  ticks=bounds, boundaries=bounds, format="%d", shrink=0.5, aspect=10)
cb.set_label("Count", fontsize="small", labelpad=2)
cb.ax.tick_params(labelsize="small")

cb.ax.set_visible(True)
# exact value annotations
for ri,row in enumerate(show):
    for ci,col in enumerate(row):
        pass
        if col < 30:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="white", ha="center", va="center", fontsize="small")
        else:
            ax.text(ci+.5, ri+.5, "{:d}".format(col), color="k", ha="center", va="center", fontsize="small")
        
#remaking labeles so it doesn't look squished
ax.set_xticks(np.arange(len(ak_pool))+.5)
lbls = np.asarray(ak_pool)
ax.set_xticklabels(["{}\nn = {}".format(ak, n) for ak, n in zip(ak_pool, primary_lob_n)], rotation=45, fontsize=7, ha="right")
# yticks
ax.set_yticks(np.arange(len(short_nuclei))+.5)

ax.set_yticklabels(["{}".format(bi) for bi in short_nuclei], fontsize="small")

plt.savefig(os.path.join(dst,"thalvtacomp_mean_counts.pdf"), bbox_inches = "tight")
