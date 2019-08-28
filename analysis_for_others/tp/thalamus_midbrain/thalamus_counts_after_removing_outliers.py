#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:18:37 2019

@author: wanglab
"""

import pickle, numpy as np, pandas as pd, matplotlib.pyplot as plt, sys, os
sys.path.append("/jukebox/wang/zahra/python/lightsheet_py3")
from tools.analysis.network_analysis import make_structure_objects

thal_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_thalamus/dataframe_no_prog_at_each_level.p"
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx"
structures = make_structure_objects(df_pth, remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

dst = "/home/wanglab/Desktop"
thal_data_pth = "/jukebox/wang/zahra/modeling/h129/thalamus/data_v2.p"
thal_data = pickle.load(open(thal_data_pth, "rb"), encoding = "latin1")

#%%
def get_counts_from_pickle(pth, sois, structures, df_pth = "/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx",
                           scale = 0.020):
    cells_regions = pickle.load(open(pth, "rb"), encoding = "latin1")
   
    #get densities for all the structures
    df = pd.read_excel(df_pth, index_col = None)
    df = df.drop(columns = ["Unnamed: 0"])
    df = df.sort_values(by = ["name"])
    
    #make new dict - for all brains
    cells_pooled_regions = {} #for raw counts
    volume_pooled_regions = {} #for density
    
    brains = list(cells_regions.keys())
    
    for brain in brains:    
        #make new dict - this is for EACH BRAIN
        c_pooled_regions = {}
        d_pooled_regions = {}
        
        for soi in sois:
            try:
                soi = [s for s in structures if s.name==soi][0]
                counts = [] #store counts in this list
                volume = [] #store volume in this list
                for k, v in cells_regions[brain].items():
                    if k == soi.name:
                        counts.append(v)
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi.name, "voxels_in_structure"].values[0])#*(scale_factor**3))
                progeny = [str(xx.name) for xx in soi.progeny]
                #now sum up progeny
                if len(progeny) > 0:
                    for progen in progeny:
                        for k, v in cells_regions[brain].items():
                            if k == progen:
                                counts.append(v)
                                #add to volume list from LUT
                        volume.append(df.loc[df.name == progen, "voxels_in_structure"].values[0])
                c_pooled_regions[soi.name] = np.sum(np.asarray(counts))
                d_pooled_regions[soi.name] = np.sum(np.asarray(volume))
            except:
                for k, v in cells_regions[brain].items():
                    if k == soi:
                        counts.append(v)                    
                #add to volume list from LUT
                volume.append(df.loc[df.name == soi, "voxels_in_structure"].values[0])
                c_pooled_regions[soi] = np.sum(np.asarray(counts))
                d_pooled_regions[soi] = np.sum(np.asarray(volume))
                        
        #add to big dict
        cells_pooled_regions[brain] = c_pooled_regions
        volume_pooled_regions[brain] = d_pooled_regions
            
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
        
    cell_counts_per_brain = np.asarray(cell_counts_per_brain)
    
    volume_per_brain = []
    i = []
    for k,v in volume_pooled_regions.items():
        dct = volume_pooled_regions[k]
        for j,l in dct.items():
            i.append(l)  
        volume_per_brain.append(i)
        #re-initialise for next
        i = []  
    volume_per_brain = np.asarray(volume_per_brain)*(scale**3)
    
    #calculate denisty
    density_per_brain = np.asarray([xx/volume_per_brain[i] for i, xx in enumerate(cell_counts_per_brain)])

    return brains, cell_counts_per_brain, density_per_brain

thal_sois = ["Parafascicular nucleus", "Posterior complex of the thalamus", "Posterior triangular thalamic nucleus",
        "Lateral posterior nucleus of the thalamus", "Lateral habenula", "Lateral dorsal nucleus of thalamus",
        "Central lateral nucleus of the thalamus", "Paraventricular nucleus of the thalamus", "Nucleus of reuniens",        
        "Mediodorsal nucleus of thalamus", "Ventral part of the lateral geniculate complex",
        "Ventral posterolateral nucleus of the thalamus", "Ventral posteromedial nucleus of the thalamus",         
        "Submedial nucleus of the thalamus", "Reticular nucleus of the thalamus",
        "Ventral medial nucleus of the thalamus", "Anteroventral nucleus of thalamus",
        "Ventral anterior-lateral complex of the thalamus"
]

thal_brains, thal_counts_per_brain, thal_density_per_brain = get_counts_from_pickle(thal_pth, thal_sois, structures)

#GET ONLY NEOCORTICAL POOLS
nc_sois = ["Thalamus"]
    
thal_brains, total_counts_per_brain, total_density_per_brain = get_counts_from_pickle(thal_pth, nc_sois, structures)

#mask to remove 3 outlier brains with  highly variable spread
curated_brains = [True, False, False, True, True, True, False, True, True, True, True, True, True, True, 
                  True, True, True, True, True, True, True, True, True]

#curate
thal_density_per_brain = thal_density_per_brain[curated_brains]
thal_brains = np.array(thal_brains)[curated_brains]
thal_counts_per_brain = thal_counts_per_brain[curated_brains]
total_counts_per_brain = total_counts_per_brain[curated_brains]
#calculate percent counts
pcounts_per_brain = (thal_counts_per_brain/total_counts_per_brain)*100

#%%

#rename nuclei
thal_sois = ["PF", "PO", "PoT", "LP","LH","LD", "CL", "PVT", "RE","MD","LGv","VPL","VPM","SMT", "RTN", "VM", "AV", "VAL"]
#first, rearrange structures in ASCENDING order (will be plotted as descending, -_-) by density and counts
density_per_brain_descending_order = np.sort(thal_density_per_brain)
order = np.argsort(np.mean(thal_density_per_brain, axis = 0))
sois_descending_density = np.array(thal_sois)[order]

pcounts_per_brain_descending_order = np.sort(pcounts_per_brain)
order = np.argsort(np.mean(pcounts_per_brain, axis = 0))
sois_descending_pcounts = np.array(thal_sois)[order]

#boxplots of densitys
plt.figure(figsize=(5,6))
plt.boxplot(density_per_brain_descending_order, vert = False, labels = sois_descending_density, sym = "", showcaps = False)
ngroup = len(density_per_brain_descending_order.T)
for i in range(ngroup):
    plt.scatter(density_per_brain_descending_order[:,i], y=np.ones(len(density_per_brain_descending_order[:,0]))*i+1, color = "k", s = 7)
plt.xlabel("Density (cells/$mm^3$)")
plt.ylabel("Thalamic nuclei")
plt.savefig(os.path.join(dst, "thal_density_boxplots.pdf"), bbox_inches = "tight")

#boxplots of percent counts
plt.figure(figsize=(5,6))
plt.boxplot(pcounts_per_brain_descending_order, vert = False, labels = sois_descending_pcounts, sym = "", showcaps = False)
ngroup = len(pcounts_per_brain_descending_order.T)
for i in range(ngroup):
    plt.scatter(pcounts_per_brain_descending_order[:,i], 
                y=np.ones(len(pcounts_per_brain_descending_order[:,0]))*i+1, color = "k", s = 10)
plt.xlabel("% of total thalamic cells")
plt.ylabel("Thalamic nuclei")
plt.savefig(os.path.join(dst, "thal_pcounts_boxplots.pdf"), bbox_inches = "tight")