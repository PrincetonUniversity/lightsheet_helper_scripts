#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:38:54 2019

@author: wanglab
"""

import pickle, numpy as np, pandas as pd, matplotlib.pyplot as plt, os
from tools.analysis.network_analysis import make_structure_objects
from scipy.stats import median_absolute_deviation as mad

#pooling regions
ann_pth = "/jukebox/LightSheetTransfer/atlas/annotation_sagittal_atlas_20um_iso.tif"
structures = make_structure_objects("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", 
                                    remove_childless_structures_not_repsented_in_ABA = True, ann_pth=ann_pth)

dst = "/home/wanglab/Desktop"
#%%
#nc
nc_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts/nc_dataframe_no_prog_at_each_level.p"
thal_pth = "/jukebox/wang/pisano/tracing_output/antero_4x_analysis/201903_antero_pooled_cell_counts_thalamus/dataframe_no_prog_at_each_level.p"

def get_counts_from_pickle(pth, sois):
    cells_regions = pickle.load(open(pth, "rb"), encoding = "latin1")
   
    #get densities for all the structures
    df = pd.read_excel("/jukebox/LightSheetTransfer/atlas/ls_id_table_w_voxelcounts.xlsx", index_col = None)
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
    return brains, cell_counts_per_brain

#GET ONLY NEOCORTICAL POOLS
nc_sois = ["Infralimbic area", "Prelimbic area", "Anterior cingulate area", "Frontal pole, cerebral cortex", "Orbital area", 
            "Gustatory areas", "Agranular insular area", "Visceral area", "Somatosensory areas", "Somatomotor areas",
            "Retrosplenial area", "Posterior parietal association areas", "Visual areas", "Temporal association areas",
            "Auditory areas", "Ectorhinal area", "Perirhinal area"]
    
nc_brains, nc_counts_per_brain = get_counts_from_pickle(nc_pth, nc_sois)

thal_sois = ["Ventral group of the dorsal thalamus", "Subparafascicular nucleus", "Subparafascicular area",
          "Peripeduncular nucleus", "Geniculate group, dorsal thalamus", "Lateral group of the dorsal thalamus",
          "Anterior group of the dorsal thalamus", "Medial group of the dorsal thalamus", "Midline group of the dorsal thalamus",
          "Intralaminar nuclei of the dorsal thalamus", "Reticular nucleus of the thalamus", "Geniculate group, ventral thalamus",
          "Epithalamus"]

thal_brains, thal_counts_per_brain = get_counts_from_pickle(thal_pth, thal_sois)
#%%
#nc
fig = plt.figure()
ax = fig.add_axes([.4,.1,.5,.8])

nc_total_counts = np.sum(nc_counts_per_brain, axis=1)
ax.hist(nc_total_counts, bins=20, color="green")
ax.set_xticks(np.arange(0,np.max(nc_total_counts), 2000)+.5)
ax.set_xticklabels(np.arange(0,np.max(nc_total_counts), 2000), rotation=30, fontsize="x-small")
ax.set_title("Total neocortical counts")
ax.set_xlabel("Neocortical cell counts")

plt.savefig(os.path.join(dst, "nc_histogram.pdf"), bbox_inches = "tight")

#thal
fig = plt.figure()
ax = fig.add_axes([.4,.1,.5,.8])

thal_total_counts = np.sum(thal_counts_per_brain, axis=1)
ax.hist(thal_total_counts, bins=20, color="green")
ax.set_xticks(np.arange(0,np.max(thal_total_counts), 200)+.5)
ax.set_xticklabels(np.arange(0,np.max(thal_total_counts), 200), rotation=30, fontsize="xx-small")
ax.set_title("Total thalamic counts")
ax.set_xlabel("Thalamic cell counts")

plt.savefig(os.path.join(dst, "thal_histogram.pdf"), bbox_inches = "tight")

#%%
df = pd.DataFrame()

est_std = 0.6745 #normal mad to get estimated standard dev
d = 3 #number of decimals
df["total_counts"] = nc_total_counts
df.index = nc_brains

df_stats = pd.DataFrame()
df_stats["mean"] = pd.Series(np.mean(nc_total_counts)).round(d)
df_stats["median"] = pd.Series(np.median(nc_total_counts)).round(d)
df_stats["std"] = pd.Series(np.std(nc_total_counts)).round(d)
df_stats["est_std"] = pd.Series(mad(nc_total_counts)/est_std).round(d)
df_stats = df_stats.T
df_stats.columns = ["total_counts"]

df_total = df.append(df_stats)
df_total.to_csv(os.path.join(dst, "nc_total_count_stats.csv"))

df = pd.DataFrame()

df["total_counts"] = thal_total_counts
df.index = thal_brains

df_stats = pd.DataFrame()
df_stats["mean"] = pd.Series(np.mean(thal_total_counts)).round(d)
df_stats["median"] = pd.Series(np.median(thal_total_counts)).round(d)
df_stats["std"] = pd.Series(np.std(thal_total_counts)).round(d)
df_stats["est_std"] = pd.Series(mad(thal_total_counts)/est_std).round(d)
df_stats = df_stats.T
df_stats.columns = ["total_counts"]

df_total = df.append(df_stats)
df_total.to_csv(os.path.join(dst, "thal_total_count_stats.csv"))